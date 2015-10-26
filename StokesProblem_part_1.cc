#include "StokesProblem.hh"
#include <cmath>
#include <algorithm>
#include <fstream>
#include <iostream>
#include <set>
#include <string>
#include <tuple>


template <int dim>
StokesProblem<dim>::StokesProblem(bool verbose, EXAMPLE example,
                                  QUADRATURE quadrature, REFINEMENT refinement, unsigned int max_degree):
  dof_handler(triangulation),
  verbose(verbose),
  example(example),
  refinement(refinement),
  max_degree(max_degree)
{
  // Set the exact solution used in the Dirichlet boundary condition and to
  // compute the error. Also set the rhs.
  switch (example)
    {
    case (example_1) :
    {
      exact_solution = new ExactSolutionEx1<dim>();
      rhs_function = new RightHandSideEx1<dim>();

      break;
    }
    case (example_2) :
    {
      exact_solution = new ExactSolutionEx2<dim>();
      AssertThrow(false, ExcMessage("Not implemented"));

      break;
    }
    case (example_3) :
    {
      exact_solution = new ExactSolutionEx3<dim>();
      rhs_function = new RightHandSideEx3<dim>();

      break;
    }
    case (example_4) :
    {
      exact_solution = new ExactSolutionEx4<dim>();
      AssertThrow(false, ExcMessage("Not implemented"));

      break;
    }
    default :
    {
      AssertThrow(false, ExcMessage("Unknow Example"));
    }
    }

  // Set the quadrature and the fe_collection
  if (quadrature==gauss_lobatto)
    {
      for (unsigned int degree=1; degree<=max_degree; ++degree)
        {
          fe_collection.push_back (FESystem<dim>(FE_Q<dim> (QGaussLobatto<1> (degree + 2)), dim,
                                                 FE_Q<dim> (QGaussLobatto<1> (degree + 1)), 1));

          quadrature_collection.push_back(QGaussLobatto<dim> (degree+3));
          face_quadrature_collection.push_back ( QGaussLobatto<dim-1> (degree+3));

          quadrature_collection_Err.push_back( QGaussLobatto<dim> (degree+4));
          face_quadrature_collection_Err.push_back ( QGaussLobatto<dim-1> (degree+4));
        }
    }
  else
    {
      for (unsigned int degree=1; degree<=max_degree; ++degree)
        {
          fe_collection.push_back (FESystem<dim>(FE_Q<dim> (degree + 1), dim,
                                                 FE_Q<dim> (degree), 1));

          quadrature_collection.push_back(QGauss<dim> (degree+2));
          face_quadrature_collection.push_back (QGauss<dim-1> (degree+2));

          quadrature_collection_Err.push_back(QGauss<dim> (degree+3));
          face_quadrature_collection_Err.push_back (QGauss<dim-1> (degree+3));
        }
    }

  fe_collection.push_back (FESystem<dim>(FE_Nothing<dim>(), dim,
                                         FE_Nothing<dim>(), 1));
  quadrature_collection.push_back(QGauss<dim>(1));
  face_quadrature_collection.push_back (QGauss<dim-1>(1));
}


template <int dim>
StokesProblem <dim>::~StokesProblem()
{
  delete exact_solution;
  delete rhs_function;
}


template <int dim>
bool StokesProblem <dim>::sort_decreasing_order (
  const std::pair<double,DoFHandler_active_cell_iterator > &i,
  const std::pair<double,DoFHandler_active_cell_iterator > &j)
{
  return ((i.first) > (j.first));
}


template <int dim>
void StokesProblem <dim>::generate_mesh()
{
  // Generate a square/cube mesh or a L-shape mesh.
  if ((example==example_3) || (example==example_3))
    {
      GridGenerator::hyper_cube(triangulation, -1, 1);
      triangulation.refine_global(2);
    }
  else
    {
      std::vector<Point<dim>> vertices(8);

      vertices[0] = Point<dim>(-1,-1);
      vertices[1] = Point<dim>(0,-1);
      vertices[2] = Point<dim>(-1,0);
      vertices[3] = Point<dim>(0,0);
      vertices[4] = Point<dim>(1,0);
      vertices[5] = Point<dim>(-1,1);
      vertices[6] = Point<dim>(0,1);
      vertices[7] = Point<dim>(1,1);

      const unsigned int n_cells(3);
      std::vector<CellData<dim> > cell(n_cells);
      cell[0].vertices[0]=0;
      cell[0].vertices[1]=1;
      cell[0].vertices[2]=2;
      cell[0].vertices[3]=3;

      cell[1].vertices[0]=2;
      cell[1].vertices[1]=3;
      cell[1].vertices[2]=5;
      cell[1].vertices[3]=6;

      cell[2].vertices[0]=3;
      cell[2].vertices[1]=4;
      cell[2].vertices[2]=6;
      cell[2].vertices[3]=7;

      triangulation.create_triangulation(vertices,cell,SubCellData());
      triangulation.refine_global(1);
    }

  // Set the degree of the FE to the lowest possible on every cell.
  DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                  end_cell = dof_handler.end();
  for (; cell!=end_cell; ++cell)
    cell->set_active_fe_index(0);
}


template <int dim>
void StokesProblem <dim>::setup_system()
{
  system_matrix.clear();

  dof_handler.distribute_dofs(fe_collection);

  std::vector<unsigned int> block_component (dim+1, 0);
  block_component[dim]=1;
  DoFRenumbering::component_wise(dof_handler, block_component);

  constraints.clear ();
  FEValuesExtractors::Vector velocities(0);
  DoFTools::make_hanging_node_constraints(dof_handler, constraints);

  VectorTools::interpolate_boundary_values(dof_handler,0,*exact_solution,constraints,fe_collection.component_mask(velocities));

  // Since with Dirichlet velocity Bdry, pressure will be defined up to a constant,
  // in order to make the solution to be unique, we need to add an additional
  // constraint for Pressure. We choose for example the first cell of triangulation
  // and do as follow:
  DoFHandler_active_cell_iterator first_cell = dof_handler.begin_active();
  std::vector<types::global_dof_index> local_dof_indices (first_cell->get_fe().dofs_per_cell);
  first_cell->get_dof_indices(local_dof_indices);
  Point<dim> pt_ref_space = first_cell->get_fe().unit_support_point(
                              first_cell->get_fe().component_to_system_index(dim,0));
  MappingQ1<dim> mapping;
  Point<dim> pt_real_space = mapping.transform_unit_to_real_cell(first_cell,pt_ref_space);
  types::global_dof_index first_pressure_dof =
    local_dof_indices[first_cell->get_fe().component_to_system_index(dim,0)];

  // component_to_system_index: "Compute the shape function for the given vector
  // component and index."
  constraints.add_line (first_pressure_dof);
  Vector<double> values(dim+1);
  exact_solution->vector_value(pt_real_space,values);
  constraints.set_inhomogeneity (first_pressure_dof,values[dim]);

  constraints.close();

  std::vector<types::global_dof_index> dofs_per_block (2);
  DoFTools::count_dofs_per_block (dof_handler, dofs_per_block, block_component);
  const unsigned int n_u = dofs_per_block[0], n_p = dofs_per_block[1];

  BlockDynamicSparsityPattern csp(2,2);
  csp.block(0,0).reinit(n_u,n_u);
  csp.block(1,0).reinit(n_p,n_u);
  csp.block(0,1).reinit(n_u,n_p);
  csp.block(1,1).reinit(n_p,n_p);
  csp.collect_sizes();
  DoFTools::make_sparsity_pattern (dof_handler, csp, constraints, false);
  sparsity_pattern.copy_from (csp);

  system_matrix.reinit(sparsity_pattern);
  solution.reinit(2);
  solution.block(0).reinit(n_u);
  solution.block(1).reinit(n_p);
  solution.collect_sizes();
  system_rhs.reinit(2);
  system_rhs.block(0).reinit (n_u);
  system_rhs.block(1).reinit (n_p);
  system_rhs.collect_sizes ();
}


template <int dim>
void StokesProblem <dim>::assemble_system ()
{
  hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection,
                                  update_values|update_quadrature_points|update_JxW_values|update_gradients);

  FullMatrix<double> local_matrix;
  Vector<double> local_rhs;
  std::vector<types::global_dof_index> local_dof_indices;

  std::vector<Vector<double>>  rhs_values;

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  // std::vector<SymmetricTensor<2,dim> > symgrad_phi_u;
  std::vector<Tensor<2,dim>> grad_phi_u;
  std::vector<double> div_phi_u;
  std::vector<Tensor<1,dim>> phi_u;
  std::vector<double> phi_p;

  DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      const unsigned int dofs_per_cell = cell->get_fe().dofs_per_cell;
      local_matrix.reinit(dofs_per_cell, dofs_per_cell);
      local_rhs.reinit(dofs_per_cell);
      local_matrix = 0;
      local_rhs = 0;

      hp_fe_values.reinit(cell);
      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
      const std::vector<double> &JxW_values = fe_values.get_JxW_values ();
      const unsigned int n_q_points = fe_values.n_quadrature_points;

      rhs_values.resize(n_q_points, Vector<double>(dim+1));
      rhs_function->vector_value_list(fe_values.get_quadrature_points(), rhs_values);

      //symgrad_phi_u.resize(dofs_per_cell);
      grad_phi_u.resize(dofs_per_cell);
      div_phi_u.resize(dofs_per_cell);
      phi_u.resize (dofs_per_cell);
      phi_p.resize(dofs_per_cell);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          for (unsigned int k=0; k<dofs_per_cell; ++k)
            {
              //symgrad_phi_u[k] = fe_values[velocities].symmetric_gradient (k, q);
              grad_phi_u[k] = fe_values[velocities].gradient (k, q);
              div_phi_u[k] = fe_values[velocities].divergence (k, q);
              phi_u[k] = fe_values[velocities].value (k, q);
              phi_p[k] = fe_values[pressure].value (k, q);
            }

          for (unsigned int i=0; i<dofs_per_cell; ++i)
            {
              for (unsigned int j=0; j<dofs_per_cell; ++j)

                //local_matrix(i,j) += (double_contract (grad_phi_u[i], grad_phi_u[j])
                //  - div_phi_u[i] * phi_p[j] - phi_p[i] * div_phi_u[j]) * JxW_values[q];
                //local_matrix(i,j) += (2 * (symgrad_phi_u[i] * symgrad_phi_u[j]) - div_phi_u[i] * phi_p[j]
                //  - phi_p[i] * div_phi_u[j])*JxW_values[q];

                local_matrix(i,j) += (double_contract<0,0,1,1> (grad_phi_u[i], grad_phi_u[j])
                                      - div_phi_u[i] * phi_p[j]- phi_p[i] * div_phi_u[j]) * JxW_values[q];


              double fe_rhs(0.);
              for (unsigned int d=0; d<dim; ++d)
                fe_rhs += phi_u[i][d]*rhs_values[q][d];
              local_rhs[i] += fe_rhs * JxW_values[q];
            }
        }

      //local system to global system
      local_dof_indices.resize (dofs_per_cell);
      cell->get_dof_indices (local_dof_indices);
      constraints.distribute_local_to_global(local_matrix, local_rhs,
                                             local_dof_indices, system_matrix, system_rhs);
    }
}


template <int dim>
void StokesProblem <dim>::solve()
{
  SparseDirectUMFPACK A_inverse;
  A_inverse.initialize(system_matrix, SparseDirectUMFPACK::AdditionalData());
  A_inverse.vmult(solution, system_rhs);
  constraints.distribute(solution);

  solution.block (1).add (-1.0 * pressure_mean_value ());
  constraints.distribute (solution);
}

template <int dim>
double StokesProblem <dim>::pressure_mean_value () const
{

  // get pressure such that satisfies mean value property:
  hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection, update_values|update_JxW_values);
  const FEValuesExtractors::Scalar pressure (dim);

  std::vector<double> values;
  double domain_mean_val_p=0;
  double measure_domain=0;
  typename hp::DoFHandler<dim>::active_cell_iterator
  cell = dof_handler.begin_active(),
  endc = dof_handler.end();
  for (; cell!=endc; ++cell)
    {
      hp_fe_values.reinit (cell);

      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
      const std::vector<double> &JxW_values = fe_values.get_JxW_values ();
      const unsigned int n_q_points = fe_values.n_quadrature_points;
      values.resize(n_q_points);
      fe_values[pressure].get_function_values(solution, values);
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          domain_mean_val_p += values[q]*JxW_values[q];
          measure_domain += JxW_values[q];
        }//q
    }//cell

  return domain_mean_val_p / measure_domain;
}


// compute_error(): Computes the Energy norm of error which is the square
// root of (L-2 norm of error in velocity gradients + L-2 norm of pressure)
template <int dim>
void StokesProblem <dim>::compute_error()
{
  error_per_cell.reinit(triangulation.n_active_cells());
  Vect_Pressure_Err.reinit(triangulation.n_active_cells());
  Vect_grad_Velocity_Err.reinit(triangulation.n_active_cells());
  Vect_Velocity_Err.reinit(triangulation.n_active_cells());

  hp::FEValues<dim> hp_fe_values(fe_collection, quadrature_collection_Err,
                                 update_values|update_quadrature_points|update_JxW_values|update_gradients|
                                 update_hessians);
  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  std::vector<double> values;
  std::vector<Tensor<2,dim>> gradients;
  std::vector<Tensor<1,dim>> velocity_values;

  std::vector<std::vector<Tensor<1,dim>>> exact_solution_gradients;
  std::vector<Vector<double>> exact_solution_values;

  DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                  endc = dof_handler.end();

  for (unsigned int cell_index=0; cell!=endc; ++cell,++cell_index)
    {
      double subtract_p = 0.;
      double grad_u_vals = 0.;
      double u_vals = 0.;
      hp_fe_values.reinit(cell);
      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
      const std::vector<double> &JxW_values = fe_values.get_JxW_values ();
      const std::vector<Point<dim> > &quadrature_points = fe_values.get_quadrature_points();
      const unsigned int n_q_points = fe_values.n_quadrature_points;

      velocity_values.resize(n_q_points);
      gradients.resize(n_q_points);
      values.resize(n_q_points);
      exact_solution_gradients.resize(n_q_points, std::vector<Tensor<1,dim>>(dim+1));
      exact_solution_values.resize(n_q_points, Vector<double> (dim+1));

      fe_values[velocities].get_function_values(solution, velocity_values);
      fe_values[velocities].get_function_gradients(solution, gradients);
      fe_values[pressure].get_function_values(solution, values);

      exact_solution->vector_gradient_list(quadrature_points, exact_solution_gradients);
      exact_solution->vector_value_list(quadrature_points, exact_solution_values);

      for (unsigned int q=0; q<n_q_points; ++q)
        {
          values[q] -= exact_solution_values[q](dim);
          subtract_p += values[q]*values[q]* JxW_values[q];

          for (unsigned int i=0; i<dim; ++i)
            {
              velocity_values[q][i] -= exact_solution_values[q](i);
              gradients[q][i] -= exact_solution_gradients[q][i];
            }

          grad_u_vals += gradients[q].norm_square() * JxW_values[q];
          u_vals += velocity_values[q].norm_square() * JxW_values[q];
        }

      error_per_cell[cell_index] = std::sqrt(subtract_p + grad_u_vals);

      Vect_Pressure_Err[cell_index] = std::sqrt(subtract_p);
      Vect_grad_Velocity_Err[cell_index] = std::sqrt(grad_u_vals);
      Vect_Velocity_Err[cell_index] = std::sqrt(u_vals);
    }

  if (verbose)
    {
      std::cout<< std::endl;
      double L2_norm_grad_velocity_Err = Vect_grad_Velocity_Err.l2_norm();
      std::cout<< "L2_norm_grad_velocity_Err: "<< L2_norm_grad_velocity_Err << std::endl;
      std::cout<< std::endl;
      double L2_norm_velocity_Err = Vect_Velocity_Err.l2_norm();
      std::cout<< "L2_norm_velocity_Err: "<< L2_norm_velocity_Err << std::endl;
      std::cout<< std::endl;
      std::cout<< std::endl;
      double L2_norm_pressure_Err = Vect_Pressure_Err.l2_norm();
      std::cout<< "L2_norm_pressure_Err: "<< L2_norm_pressure_Err << std::endl;
      std::cout<< std::endl;
      std::cout<< std::endl;
      double L2_norm_total_Err = sqrt(std::pow(L2_norm_grad_velocity_Err,2) +
                                      std::pow (L2_norm_pressure_Err,2));
      std::cout<< "L2_norm of Total_ERROR is: "<< L2_norm_total_Err << std::endl;
      std::cout<< std::endl;
    }
}

// estimate_error(): Computes the residual based error estimator in
// hp- adaptive FEM approach
template <int dim>
void StokesProblem<dim>::estimate_error()
{
  est_per_cell.reinit(triangulation.n_active_cells());
  hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection,
                                  update_values|update_quadrature_points|update_JxW_values|update_gradients|
                                  update_hessians);
  hp::FEFaceValues<dim> hp_fe_face_values(fe_collection, face_quadrature_collection,
                                          update_JxW_values|update_gradients|update_normal_vectors);
  hp::FEFaceValues<dim> hp_neighbor_face_values(fe_collection,
                                                face_quadrature_collection, update_gradients);
  hp::FESubfaceValues<dim> hp_subface_values(fe_collection, face_quadrature_collection,
                                             update_JxW_values|update_gradients|update_normal_vectors);
  hp::FESubfaceValues<dim> hp_neighbor_subface_values(fe_collection,
                                                      face_quadrature_collection, update_gradients);

  std::vector<Tensor<1,dim> > gradients_p;
  std::vector<double> divergences;
  std::vector<Tensor<1,dim> > laplacians;

  std::vector<Tensor<2,dim> > gradients;
  std::vector<Tensor<2,dim> > neighbor_gradients;

  const FEValuesExtractors::Vector velocities (0);
  const FEValuesExtractors::Scalar pressure (dim);

  std::vector<Vector<double>> rhs_values;

  Vector<double> res_est_per_cell(triangulation.n_active_cells());
  Vector<double> Jump_est_per_cell(triangulation.n_active_cells());

  DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                  endc = dof_handler.end();
  unsigned int cell_index=0;
  for (; cell!=endc; ++cell,++cell_index)
    {
      hp_fe_values.reinit (cell);
      const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
      const std::vector<double> &JxW_values = fe_values.get_JxW_values ();
      const unsigned int n_q_points = fe_values.n_quadrature_points;

      rhs_values.resize(n_q_points, Vector<double>(dim+1));
      rhs_function->vector_value_list(fe_values.get_quadrature_points(), rhs_values);

      divergences.resize(n_q_points);
      gradients_p.resize(n_q_points) ;
      laplacians.resize(n_q_points);

      fe_values[pressure].get_function_gradients(solution, gradients_p);
      fe_values[velocities].get_function_divergences(solution, divergences);
      fe_values[velocities].get_function_laplacians(solution, laplacians);


      double term2=0;//divergence term in estimator definition
      double term1=0;// for the residual term in estimator definition
      for (unsigned int q=0; q<n_q_points; ++q)
        {
          term2 += (divergences[q])*(divergences[q])*JxW_values[q];

          for (unsigned int i=0; i<2; ++i)
            gradients_p[q][i] -= (rhs_values[q](i)+laplacians[q][i]);

          term1 += contract<0,0>(gradients_p[q],gradients_p[q])*JxW_values[q];
        }
      res_est_per_cell[cell_index] = pow((cell->diameter())/(cell->get_fe().degree), 2.0) * (term1) + term2;

      //  compute jump_est_per_cell
      double term3=0;//jump part of the estimator
      for (unsigned int face_number=0; face_number<GeometryInfo<dim>::faces_per_cell; ++face_number)
        {
          // Check that the face is not at the boundary
          if (cell->face(face_number)->at_boundary()==false)
            {
              // If the neighbor cell is active and on the same level
              if ((cell->neighbor(face_number)->active()) &&
                  (cell->neighbor(face_number)->level() == cell->level()))
                {
                  const unsigned int q_index = std::max (cell->active_fe_index(),
                                                         cell->neighbor(face_number)->active_fe_index());

                  hp_fe_face_values.reinit(cell, face_number, q_index);
                  hp_neighbor_face_values.reinit (cell->neighbor(face_number),
                                                  cell->neighbor_of_neighbor(face_number), q_index);

                  const FEFaceValues<dim> &neighbor_face_values = hp_neighbor_face_values.get_present_fe_values();
                  const FEFaceValues<dim> &fe_face_values = hp_fe_face_values.get_present_fe_values ();
                  const std::vector<double> &JxW_values = fe_face_values.get_JxW_values ();
                  const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

                  gradients.resize(n_face_q_points);
                  neighbor_gradients.resize(n_face_q_points);

                  neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
                  fe_face_values[velocities].get_function_gradients(solution, gradients);

                  std::vector<Tensor<1,dim> > jump_per_face;
                  jump_per_face.resize(n_face_q_points);
                  double jump_val=0;

                  for (unsigned int q=0; q<n_face_q_points; ++q)
                    {
                      for (unsigned int i=0; i<2; ++i)
                        for (unsigned int j=0; j<2; ++j)
                          jump_per_face[q][i] = (gradients[q][i][j]-neighbor_gradients[q][i][j]) *
                                                (fe_face_values.normal_vector(q)[j]);
                      jump_val += contract<0,0>(jump_per_face[q],jump_per_face[q])*JxW_values[q];
                    }

                  unsigned int min_fe_degree = std::min(cell->get_fe().degree,
                                                        cell->neighbor(face_number)->get_fe().degree);


                  term3 += (cell->face(face_number)->diameter())/(2.0*min_fe_degree)*jump_val;
                }
              // If the neighbor has children
              else if (cell->neighbor(face_number)->has_children() == true)
                {
                  for (unsigned int subface=0; subface< cell->face(face_number)->n_children(); ++subface)
                    {
                      const unsigned int q_index = std::max(cell->neighbor_child_on_subface(
                                                              face_number, subface)->active_fe_index(), cell->active_fe_index());


                      hp_neighbor_face_values.reinit(cell->neighbor_child_on_subface(
                                                       face_number, subface), cell->neighbor_of_neighbor(face_number), q_index);
                      hp_subface_values.reinit (cell,face_number, subface, q_index);

                      const FEFaceValues<dim> &neighbor_face_values = hp_neighbor_face_values.get_present_fe_values();
                      const FESubfaceValues<dim> &fe_subface_values = hp_subface_values.get_present_fe_values();
                      const std::vector<double> &JxW_values = fe_subface_values.get_JxW_values ();
                      const unsigned int n_subface_q_points = fe_subface_values.n_quadrature_points;

                      gradients.resize(n_subface_q_points);
                      neighbor_gradients.resize(n_subface_q_points);

                      neighbor_face_values[velocities].get_function_gradients(solution, neighbor_gradients);
                      fe_subface_values[velocities].get_function_gradients(solution, gradients);

                      std::vector<Tensor<1,dim>> jump_per_subface;
                      jump_per_subface.resize(n_subface_q_points);

                      double jump_val=0;
                      for (unsigned int q=0; q<n_subface_q_points; ++q)
                        {
                          for (unsigned int i=0; i<2; ++i)
                            for (unsigned int j=0; j<2; ++j)
                              jump_per_subface[q][j] += (gradients[q][i][j] - neighbor_gradients[q][i][j]) *
                                                        (fe_subface_values.normal_vector(q)[j]);
                          jump_val += contract<0,0>(jump_per_subface[q],jump_per_subface[q])*(JxW_values[q]);
                        }

                      unsigned int min_fe_degree = std::min(cell->get_fe().degree,
                                                            cell->neighbor_child_on_subface(face_number,subface)->get_fe().degree);

                      term3 +=(cell->face(face_number)->child(subface)->diameter()) / (2.0 *
                              min_fe_degree)*jump_val;
                    }
                }
              // The neighbor is coarser
              else
                {
                  Assert(cell->neighbor_is_coarser(face_number),
                         ExcMessage("Problem while computing the error estimator."));
                  const unsigned int q_index = std::max(cell->active_fe_index(),
                                                        cell->neighbor(face_number)->active_fe_index());
                  hp_fe_face_values.reinit(cell, face_number,q_index);
                  hp_neighbor_subface_values.reinit(cell->neighbor(face_number),
                                                    cell->neighbor_of_coarser_neighbor(face_number).first,
                                                    cell->neighbor_of_coarser_neighbor(face_number).second,q_index);

                  const FEFaceValues<dim> &fe_face_values = hp_fe_face_values.get_present_fe_values();
                  const FESubfaceValues<dim> &neighbor_subface_values =
                    hp_neighbor_subface_values.get_present_fe_values();

                  const std::vector<double> &JxW_values = fe_face_values.get_JxW_values();

                  const unsigned int n_face_q_points = fe_face_values.n_quadrature_points;

                  gradients.resize(n_face_q_points);
                  neighbor_gradients.resize(n_face_q_points);

                  neighbor_subface_values[velocities].get_function_gradients(solution, neighbor_gradients);
                  fe_face_values[velocities].get_function_gradients(solution, gradients);

                  std::vector<Tensor<1,dim> > jump_per_face;
                  jump_per_face.resize(n_face_q_points);
                  double jump_val=0;
                  for (unsigned int q=0; q<n_face_q_points; ++q)
                    {
                      for (unsigned int i=0; i<2; ++i)
                        for (unsigned int j=0; j<2; ++j)
                          jump_per_face[q][i] += (gradients[q][i][j] - neighbor_gradients[q][i][j]) *
                                                 (fe_face_values.normal_vector(q)[j]);
                      jump_val += contract<0,0>(jump_per_face[q],jump_per_face[q])*JxW_values[q];
                    }


                  unsigned int min_fe_degree = std::min(cell->get_fe().degree,
                                                        cell->neighbor(face_number)->get_fe().degree);

                  term3 += (cell->face(face_number)->diameter())/(2.0*min_fe_degree)*jump_val;
                }
            }

        }
      Jump_est_per_cell(cell_index) = term3;
      est_per_cell(cell_index) = sqrt(Jump_est_per_cell(cell_index)+res_est_per_cell(cell_index));
    }
}

// get_patch_around_cell() : For each cell returns a vector of cells
// which are located around that cell
template <int dim>
std::vector<typename StokesProblem<dim>::DoFHandler_active_cell_iterator>
StokesProblem<dim>::get_patch_around_cell(const DoFHandler_active_cell_iterator &cell)
{
  std::vector<DoFHandler_active_cell_iterator> patch;
  std::set<DoFHandler_active_cell_iterator> cells_done;

  patch.push_back (cell);
  cells_done.insert(cell);
//  i counter for the number of patch layers ... n_layers=1 here (1 level of patch around cell)
  for (unsigned int i=0; i<1; ++i)
    {
      const unsigned int patch_size = patch.size();
      for (unsigned int j=0; j<patch_size; ++j)
        {
          for (unsigned int face_number=0; face_number< GeometryInfo<dim>::faces_per_cell; ++face_number)
            {
              if (patch[j]->face(face_number)->at_boundary()==false)
                {
                  if (patch[j]->face(face_number)->has_children() == false)
                    {
                      DoFHandler_active_cell_iterator cell = patch[j]->neighbor(face_number);
                      if (cells_done.count(cell)==0)
                        {
                          patch.push_back(cell);
                          cells_done.insert(cell);
                        }
                    }
                  else
                    {
                      for (unsigned int subface=0; subface< patch[j]->face(face_number)->n_children(); ++subface)
                        {
                          DoFHandler_active_cell_iterator child_cell =
                            patch[j]->neighbor_child_on_subface (face_number, subface);
                          if (cells_done.count(child_cell)==0)
                            {
                              patch.push_back(child_cell);
                              cells_done.insert(child_cell);
                            }
                        }
                    }
                }
            }
        }
    }

  return patch;
}

// In order to build triangulation from patch cells, first we need to
// return the vector of the coarsest common level cells as follows
template <int dim>
std::vector<typename StokesProblem<dim>::DoFHandler_cell_iterator>
StokesProblem<dim>::get_cells_at_coarsest_common_level (
  const std::vector<DoFHandler_active_cell_iterator> &patch)
{
  Assert (patch.size() > 0, ExcMessage("vector containing patch cells should not be an empty vector!"));
  unsigned int min_level = static_cast<unsigned int> (patch[0]->level());
  unsigned int max_level = static_cast<unsigned int> (patch[0]->level());
  for (unsigned int i=0; i<patch.size(); ++i)
    {
      min_level = std::min (min_level, static_cast<unsigned int> (patch[i]->level()) );
      max_level = std::max (max_level, static_cast<unsigned int> (patch[i]->level()) );
    }

  std::set<DoFHandler_cell_iterator>  uniform_cells;

  typename std::vector<DoFHandler_active_cell_iterator>::const_iterator  patch_cell;

  for (patch_cell=patch.begin(); patch_cell!=patch.end () ; ++patch_cell)
    {
      if (static_cast<unsigned int>((*patch_cell)->level()) == min_level)
        uniform_cells.insert (*patch_cell);
      else
        {
          DoFHandler_cell_iterator parent = *patch_cell;

          while (static_cast<unsigned int> (parent->level()) > min_level)
            parent = parent-> parent();
          uniform_cells.insert (parent);
        }
    }

  return std::vector<DoFHandler_cell_iterator> (uniform_cells.begin(),
                                                uniform_cells.end());
}

// To solve local variational problems on each patch, we have to build the corresponding triangulation on
// each patch
template <int dim>
void StokesProblem<dim>::build_triangulation_from_patch(
  const std::vector<DoFHandler_active_cell_iterator> &patch,
  Triangulation<dim> &local_triangulation, unsigned int &level_h_refine,
  unsigned int &level_p_refine, std::map<Triangulation_active_cell_iterator,
  DoFHandler_active_cell_iterator> &patch_to_global_tria_map)
{
  std::vector<DoFHandler_cell_iterator> uniform_cells =
    get_cells_at_coarsest_common_level (patch);

  level_h_refine=static_cast<unsigned int> (patch[0]->level());
  level_p_refine=static_cast<unsigned int> (patch[0]->active_fe_index());

  local_triangulation.clear();
  std::vector<Point<dim> > vertices;
  const unsigned int n_uniform_cells=uniform_cells.size();
  std::vector<CellData<dim> > cells(n_uniform_cells);
  unsigned int k=0;// for enumerating cells
  unsigned int i=0;// for enumerating vertices

  typename std::vector<DoFHandler_cell_iterator>::const_iterator uniform_c;

  for (uniform_c=uniform_cells.begin(); uniform_c!=uniform_cells.end(); ++uniform_c)
    {
      bool repeat_vertex;
      for (unsigned int j=0;  j< GeometryInfo<dim>::vertices_per_cell; ++j)
        {
          Point<dim> position=(*uniform_c)->vertex (j);
          repeat_vertex=false;

          for (unsigned int m=0; m<i; ++m)
            {
              if (position == vertices[m])
                {
                  repeat_vertex=true;
                  cells[k].vertices[j]=m ;
                  break;
                }
            }
          if (repeat_vertex==false)
            {
              vertices.push_back(position);
              cells[k].vertices[j]=i;
              i=i+1;
            }

        }//for vertices_per_cell
      k=k+1;
    }

  local_triangulation.create_triangulation(vertices,cells,SubCellData());
  Assert (local_triangulation.n_active_cells() == uniform_cells.size(), ExcInternalError());

  local_triangulation.clear_user_flags ();
  unsigned int index=0;
  std::map<Triangulation_cell_iterator, DoFHandler_cell_iterator> patch_to_global_tria_map_tmp;
  for (Triangulation_cell_iterator coarse_cell_t = local_triangulation.begin();
       coarse_cell_t != local_triangulation.end(); ++coarse_cell_t, ++index)
    {
      patch_to_global_tria_map_tmp.insert (std::make_pair(coarse_cell_t, uniform_cells[index]));
      AssertThrow ((std::fabs(coarse_cell_t->center()(0) - uniform_cells[index]->center()(0))<1e-16 &&
                    std::fabs(coarse_cell_t->center()(1) - uniform_cells[index]->center()(1)) <1e-16),
                   ExcInternalError());
    }

  bool refinement_necessary;
  do
    {
      refinement_necessary = false;
      for (Triangulation_active_cell_iterator cell_tt = local_triangulation.begin_active();
           cell_tt != local_triangulation.end(); ++cell_tt)
        {
          if (patch_to_global_tria_map_tmp[cell_tt]->has_children())
            {
              cell_tt -> set_refine_flag();
              refinement_necessary = true;
            }
          else for (unsigned int i=0; i<patch.size(); ++i)
              {
                if (patch_to_global_tria_map_tmp[cell_tt]==patch[i])
                  {
                    cell_tt->set_user_flag();
                    break;
                  }
              }
        }

      if (refinement_necessary)
        {
          local_triangulation.execute_coarsening_and_refinement ();

          for (Triangulation_cell_iterator cell_ttt = local_triangulation.begin();
               cell_ttt != local_triangulation.end(); ++cell_ttt)
            {

              if (patch_to_global_tria_map_tmp.find(cell_ttt)!=patch_to_global_tria_map_tmp.end())
                {
                  if (cell_ttt-> has_children())
                    {
                      // Note: Since the cell got children, then it should not be in the map anymore...
                      // children may be added into the map, instead

                      // these children may not yet be in the map
                      for (unsigned int c=0; c< cell_ttt ->n_children(); ++c)
                        {
                          if (patch_to_global_tria_map_tmp.find(cell_ttt->child(c)) ==
                              patch_to_global_tria_map_tmp.end())
                            {
                              patch_to_global_tria_map_tmp.insert (std::make_pair(cell_ttt ->child(c),
                                                                                  patch_to_global_tria_map_tmp[cell_ttt]->child(c)));

                              AssertThrow( (std::fabs (cell_ttt ->child(c)->center()(0) -
                                                       patch_to_global_tria_map_tmp[cell_ttt]->child(c)->center()(0)) < 1e-16 &&
                                            std::fabs (cell_ttt ->child(c)->center()(1) -
                                                       patch_to_global_tria_map_tmp[cell_ttt]->child(c)->center()(1)) < 1e-16),
                                           ExcInternalError());
                            }
                        }
                      patch_to_global_tria_map_tmp.erase(cell_ttt);
                    }
                }
            }
        }

    }
  while (refinement_necessary);

  typename std::map<Triangulation_cell_iterator,DoFHandler_cell_iterator>::iterator map_tmp_it =
    patch_to_global_tria_map_tmp.begin(),map_tmp_end = patch_to_global_tria_map_tmp.end();

  for (; map_tmp_it!=map_tmp_end; ++map_tmp_it)
    patch_to_global_tria_map[map_tmp_it->first] = map_tmp_it->second;
}


// mark cells in the "local_triangulation" that exist in the patch: go
// through all cells in the "local_triangulation" and see whether the
// corresponding cell in the global triangulation is part of
// the 'patch' list of cells
template <int dim>
void StokesProblem<dim>::set_active_fe_indices(hp::DoFHandler<dim> &local_dof_handler,
                                               std::map<Triangulation_active_cell_iterator, DoFHandler_active_cell_iterator>
                                               &patch_to_global_tria_map)
{
  DoFHandler_active_cell_iterator patch_cell = local_dof_handler.begin_active(),
                                  end_patch_cell = local_dof_handler.end();
  for (; patch_cell!=end_patch_cell; ++patch_cell)
    {
      if (patch_cell->user_flag_set()==true)
        {
          DoFHandler_active_cell_iterator global_cell = patch_to_global_tria_map[patch_cell];

          patch_cell->set_active_fe_index(global_cell->active_fe_index());
        }
      else if (patch_cell->user_flag_set()==false)
        {
          // which assigns FE_Nothing for the cells out of patch
          patch_cell->set_active_fe_index (max_degree);
        }
      else
        Assert (false, ExcNotImplemented());
    }
}


