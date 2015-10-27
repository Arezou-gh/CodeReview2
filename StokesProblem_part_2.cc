
template <int dim>
void StokesProblem<dim>::patch_output (unsigned int patch_number,
                                       const unsigned int cycle, hp::DoFHandler<dim> &local_dof_handler,
                                       BlockVector<double> &local_solu)
{
  std::vector<std::string> solution_names;
  solution_names.push_back("patch_x_velocity");
  solution_names.push_back("patch_y_velocity");
  solution_names.push_back("patch_pressure");

  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation;
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  DataOut<dim,hp::DoFHandler<dim> > patch_data_out;
  patch_data_out.attach_dof_handler (local_dof_handler);


  patch_data_out.add_data_vector(local_solu, solution_names,
                                 DataOut<dim,hp::DoFHandler<dim> >::type_dof_data, data_component_interpretation);
  patch_data_out.build_patches();

  std::string filename = "patch_solution-" +
                         Utilities::int_to_string (cycle, 2) +
                         +"-"+Utilities::int_to_string (patch_number, 2) +".vtu";
  std::ofstream output (filename.c_str());
  patch_data_out.write_vtu (output);
}


template <int dim>
void StokesProblem<dim>::p_refinement(hp::DoFHandler<dim> &local_dof_handler,
                                      std::map<Triangulation_active_cell_iterator, DoFHandler_active_cell_iterator>
                                      &patch_to_global_tria_map, unsigned int level_p_refine, BlockVector<double> &local_solu)
{
  bool need_to_refine = false;

  DoFHandler_active_cell_iterator act_patch_cell = local_dof_handler.begin_active(),
                                  act_end_patch_cell = local_dof_handler.end();
  for (; act_patch_cell!=act_end_patch_cell; ++act_patch_cell)
    {
      DoFHandler_active_cell_iterator global_cell = patch_to_global_tria_map[act_patch_cell];

      // Refine the cell if the refinement level is below the refinement of the
      // ``main'' cell + 1 and below the maximum level of refinement.
      if ((global_cell->active_fe_index()+1) < (fe_collection.size()-1) &&
          global_cell->active_fe_index() < (level_p_refine+1))
        need_to_refine = true;
    }

  if (need_to_refine==true)
    {
      //SolutionTransfer<dim,Vector<double>, hp::DoFHandler<dim> > aaa(local_dof_handler);
      SolutionTransfer<dim,BlockVector<double>, hp::DoFHandler<dim> > solution_transfer(local_dof_handler);
      solution_transfer.prepare_for_pure_refinement();

      DoFHandler_active_cell_iterator active_patch_cell = local_dof_handler.begin_active();
      for ( ; active_patch_cell!=local_dof_handler.end(); ++active_patch_cell)
        {
          // since here the fe_collection.size()=7 (i.e., 6 indices in total), and
          // we also know the last index hold for the fe_nothing FE, therefore we will
          // set_active_fe_index for the cells up to index 4. (the reason is that for
          // example we cannot exceed the polynomial degree for the index 5...it is
          // already the last index before the fe_nothing (index=6))
          // It is also worth to see how they did p-refinement in step-27.
          // if  (cell->active_fe_index()+1 < fe_collection.size()))
          DoFHandler_active_cell_iterator global_cell = patch_to_global_tria_map[active_patch_cell];
          active_patch_cell->set_active_fe_index(global_cell->active_fe_index()+1);
        }
      set_active_fe_indices(local_dof_handler,patch_to_global_tria_map);
      local_dof_handler.distribute_dofs (fe_collection);

      std::vector<unsigned int> block_component_patch(dim+1, 0);
      block_component_patch[dim] = 1;
      DoFRenumbering::component_wise(local_dof_handler, block_component_patch);

      std::vector<types::global_dof_index> dofs_per_block_patch (2);
      DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
      // resize the vector temp to the correct size
      BlockVector<double> temp (dofs_per_block_patch);
      solution_transfer.refine_interpolate(local_solu , temp);
      local_solu = temp;
    }
}


template <int dim>
void StokesProblem<dim>::h_refinement(Triangulation<dim> &local_triangulation,
                                      hp::DoFHandler<dim> &local_dof_handler, unsigned int level_h_refine,
                                      BlockVector<double> &local_solu)
{
  bool need_to_refine = false;
  DoFHandler_active_cell_iterator act_patch_cell = local_dof_handler.begin_active(),
                                  act_end_patch_cell = local_dof_handler.end();
  for (; act_patch_cell!=act_end_patch_cell; ++act_patch_cell)
    {
      if (static_cast<unsigned int> (act_patch_cell->level()) < (level_h_refine+1))
        {
          need_to_refine = true;
          act_patch_cell->set_refine_flag();
        }
    }

  if (need_to_refine==true)
    {
      // user flags will be overwritten by
      // execute_coarsening_and_refinement. save their values into
      // the material_id, since that one not only survives
      // refinement but is also inherited to the children
      for (DoFHandler_active_cell_iterator patch_cell = local_dof_handler.begin_active();
           patch_cell != local_dof_handler.end(); ++patch_cell)
        {
          if (patch_cell->user_flag_set())
            patch_cell->set_material_id (1);
          else
            patch_cell->set_material_id (0);
        }

      local_triangulation.prepare_coarsening_and_refinement ();
      SolutionTransfer<dim,BlockVector<double>, hp::DoFHandler<dim>> solution_transfer(local_dof_handler);
      solution_transfer.prepare_for_pure_refinement();

      local_triangulation.execute_coarsening_and_refinement ();

      // get user flags back out of the material_id field
      for (DoFHandler_cell_iterator patch_cell = local_dof_handler.begin ();
           patch_cell!=local_dof_handler.end(); ++patch_cell)
        {
          if (patch_cell->material_id() == 1)
            patch_cell->set_user_flag();
          else
            patch_cell->clear_user_flag();
        }
      local_dof_handler.distribute_dofs(fe_collection);

      std::vector<unsigned int> block_component_patch(dim+1, 0);
      block_component_patch[dim] = 1;
      DoFRenumbering::component_wise(local_dof_handler, block_component_patch);

      std::vector<types::global_dof_index> dofs_per_block_patch(2);
      DoFTools::count_dofs_per_block (local_dof_handler, dofs_per_block_patch, block_component_patch);
      // resize the vector temp to the correct size
      BlockVector<double> tmp (dofs_per_block_patch);
      solution_transfer.refine_interpolate(local_solu, tmp);
      local_solu = tmp;
    }
}


template <int dim>
void StokesProblem<dim>::patch_assemble_system(hp::DoFHandler<dim> const &local_dof_handler,
                                               ConstraintMatrix const &constraints_patch, BlockVector<double> const &local_solu,
                                               BlockSparseMatrix<double> &patch_system, BlockVector<double> &patch_rhs)
{
  hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection,
                                  update_values|update_quadrature_points|update_JxW_values|update_gradients|
                                  update_hessians);

  FullMatrix<double> local_matrix_patch;
  Vector<double> local_rhs_patch;
  std::vector<types::global_dof_index> local_dof_indices;

  std::vector<Vector<double>> rhs_values;

  FEValuesExtractors::Vector velocities (0);
  FEValuesExtractors::Scalar pressure (dim);

  std::vector<Tensor<2,dim>> grad_phi_u;
  std::vector<double> div_phi_u;
  std::vector<Tensor<1,dim>> phi_u;
  std::vector<double> phi_p;

  std::vector<Tensor<1,dim>> gradients_p;
  std::vector<double> divergences;
  std::vector<Tensor<1,dim>> laplacians;

  std::vector<double> values;
  std::vector<Tensor<2,dim>> gradients;

  DoFHandler_active_cell_iterator act_patch_cell = local_dof_handler.begin_active(),
                                  act_end_patch_cell = local_dof_handler.end();
  for (; act_patch_cell!=act_end_patch_cell; ++act_patch_cell)
    {
      const unsigned int dofs_per_cell = act_patch_cell->get_fe().dofs_per_cell;
      if (dofs_per_cell!=0)
        {
          local_matrix_patch.reinit (dofs_per_cell, dofs_per_cell);
          local_rhs_patch.reinit (dofs_per_cell);

          hp_fe_values.reinit (act_patch_cell);
          const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
          const std::vector<double> &JxW_values = fe_values.get_JxW_values ();
          const unsigned int n_q_points = fe_values.n_quadrature_points;

          rhs_values.resize(n_q_points, Vector<double>(dim+1));
          rhs_function->vector_value_list(fe_values.get_quadrature_points(), rhs_values);

          grad_phi_u.resize(dofs_per_cell);
          div_phi_u.resize(dofs_per_cell);
          phi_u.resize (dofs_per_cell);
          phi_p.resize(dofs_per_cell);

          divergences.resize(n_q_points);
          gradients_p.resize(n_q_points);
          laplacians.resize(n_q_points);

          fe_values[pressure].get_function_gradients(local_solu, gradients_p);
          fe_values[velocities].get_function_divergences(local_solu, divergences);
          fe_values[velocities].get_function_laplacians(local_solu, laplacians);


          for (unsigned int q=0; q<n_q_points; ++q)
            {
              Vector<double> local_rhs1;
              Vector<double> local_rhs2;
              local_rhs1.reinit (dofs_per_cell);
              local_rhs2.reinit (dofs_per_cell);
              // Optimizations
              const double JxW_val(JxW_values[q]);
              const double div(divergences[q]);
              std::vector<double> alpha(dim,0.);
              for (unsigned int d=0; d<dim; ++d)
                alpha[d] = rhs_values[q][d] + laplacians[q][d]-gradients_p[q][d];

              for (unsigned int k=0; k<dofs_per_cell; ++k)
                {
                  grad_phi_u[k] = fe_values[velocities].gradient (k, q);
                  phi_u[k] = fe_values[velocities].value (k, q);
                  phi_p[k] = fe_values[pressure].value (k, q);
                }


              for (unsigned int i=0; i<dofs_per_cell; ++i)
                {
                  for (unsigned int j=0; j<dofs_per_cell; ++j)
                    local_matrix_patch(i,j) +=
                      (double_contract<0,0,1,1>(grad_phi_u[i],grad_phi_u[j]) +
                       (phi_p[i]*phi_p[j])) * JxW_val;

                  for (unsigned int d=0; d<dim; ++d)
                    local_rhs1[i] += alpha[d]*phi_u[i][d];
                  local_rhs1[i] *= JxW_val;

                  local_rhs2[i] = phi_p[i]*div*JxW_val;
                  local_rhs_patch[i] += local_rhs1(i)+local_rhs2(i);
                }
            }

          local_dof_indices.resize (dofs_per_cell);
          act_patch_cell->get_dof_indices (local_dof_indices);
          constraints_patch.distribute_local_to_global(local_matrix_patch, local_rhs_patch,
                                                       local_dof_indices, patch_system, patch_rhs);
        }
    }
}


template <int dim>
void StokesProblem<dim>::patch_solve(hp::DoFHandler<dim> &local_dof_handler,
                                     unsigned int patch_number, unsigned int cycle, BlockVector<double> &local_solu,
                                     double &conv_est, double &workload_num)
{
  // setup_patch_system and patch_rhs
  workload_num = local_dof_handler.n_dofs();
  if (verbose==true)
    patch_output(patch_number, cycle, local_dof_handler, local_solu);
  ConstraintMatrix constraints_patch;
  FEValuesExtractors::Vector velocities(0);
  FEValuesExtractors::Scalar pressure (dim);
  DoFTools::make_hanging_node_constraints(local_dof_handler, constraints_patch);

  // Zero_Bdry_Condition_on_Patch
  DoFHandler_active_cell_iterator act_patch_cell = local_dof_handler.begin_active(),
                                  act_end_patch_cell = local_dof_handler.end();
  for (; act_patch_cell!=act_end_patch_cell; ++act_patch_cell)
    {
      std::vector<types::global_dof_index> local_face_dof_indices((act_patch_cell->get_fe()).dofs_per_face);
      if (act_patch_cell->user_flag_set() == true)
        {
          for (unsigned int f=0; f<GeometryInfo<dim>::faces_per_cell; ++f)
            {
              bool face_is_on_patch_Bdry = false;
              if (act_patch_cell->face(f)->at_boundary())
                face_is_on_patch_Bdry = true;
              else
                {
                  if (act_patch_cell->neighbor(f)->has_children() == true)
                    {
                      for (unsigned int sf=0; sf<act_patch_cell->face(f)->n_children(); ++sf)
                        if (act_patch_cell->neighbor_child_on_subface(f, sf)->user_flag_set() == false)
                          {
                            face_is_on_patch_Bdry = true;
                            break;
                          }
                    }
                  else
                    {
                      if (act_patch_cell->neighbor(f)->user_flag_set() == false)
                        face_is_on_patch_Bdry = true;
                    }
                }

              if (face_is_on_patch_Bdry)
                {
                  act_patch_cell->face(f)->get_dof_indices (local_face_dof_indices,
                                                            act_patch_cell->active_fe_index());
                  for (unsigned int i=0; i<local_face_dof_indices.size(); ++i)
                    // system_to_component_index: "Compute vector component and index
                    // of this shape function within the shape functions corresponding
                    // to this component from the index of a shape function within this
                    // finite element"
                    if ((act_patch_cell->get_fe()).face_system_to_component_index(i).first < dim)
                      constraints_patch.add_line (local_face_dof_indices[i]);
                }
            }
        }
    }
  constraints_patch.close();


  std::vector<unsigned int> block_component_patch (dim+1, 0);
  block_component_patch[dim]=1;
  DoFRenumbering::component_wise(local_dof_handler, block_component_patch);
  std::vector<types::global_dof_index> dofs_per_block_patch (2);
  DoFTools::count_dofs_per_block(local_dof_handler, dofs_per_block_patch, block_component_patch);

  BlockDynamicSparsityPattern csp (dofs_per_block_patch, dofs_per_block_patch);
  BlockSparsityPattern sparsity_pattern_patch;
  DoFTools::make_sparsity_pattern (local_dof_handler, csp, constraints_patch, false);
  sparsity_pattern_patch.copy_from(csp);

  BlockSparseMatrix<double> patch_system (sparsity_pattern_patch);
  BlockVector<double> patch_solution (dofs_per_block_patch);
  BlockVector<double> patch_rhs (dofs_per_block_patch);

  // assemble patch_system and patch_rhs
  patch_assemble_system(local_dof_handler, constraints_patch, local_solu, patch_system,
                        patch_rhs);

  // iterative solver
  double tolerance = 1e-12;
  SolverControl solver_control_stiffness (patch_rhs.block(0).size(),
                                          tolerance*patch_rhs.block(0).l2_norm());
  SolverCG<> cg_stiff (solver_control_stiffness);

  PreconditionSSOR<> preconditioner_stiffness;
  preconditioner_stiffness.initialize(patch_system.block(0,0), 1.2);
  cg_stiff.solve (patch_system.block(0,0), patch_solution.block(0), patch_rhs.block(0),
                  preconditioner_stiffness);

  SolverControl solver_control_mass (patch_rhs.block(1).size(),
                                     tolerance*patch_rhs.block(1).l2_norm());
  SolverCG<> cg_mass (solver_control_mass);

  PreconditionSSOR<> preconditioner_mass;
  preconditioner_mass.initialize(patch_system.block(1,1), 1.2);
  cg_mass.solve (patch_system.block(1,1), patch_solution.block(1), patch_rhs.block(1),
                 preconditioner_mass);
  constraints_patch.distribute(patch_solution);

  // get the L2 norm of the gradient of velocity solution and pressure value
  double pressure_val=0;
  double grad_u_val=0;
  double solu_norm_per_patch=0.;
  std::vector<double> values;
  std::vector<Tensor<2,dim>> gradients;
  hp::FEValues<dim> hp_fe_values (fe_collection, quadrature_collection,
                                  update_values|update_quadrature_points|update_JxW_values|update_gradients|
                                  update_hessians);

  act_patch_cell = local_dof_handler.begin_active();
  for (; act_patch_cell!=act_end_patch_cell; ++act_patch_cell)
    {
      const unsigned int dofs_per_cel = act_patch_cell->get_fe().dofs_per_cell;
      if (dofs_per_cel!=0)
        {
          hp_fe_values.reinit (act_patch_cell);
          const FEValues<dim> &fe_values = hp_fe_values.get_present_fe_values ();
          const std::vector<double> &JxW_values = fe_values.get_JxW_values ();
          const unsigned int n_q_points = fe_values.n_quadrature_points;

          gradients.resize(n_q_points);
          values.resize(n_q_points);

          fe_values[velocities].get_function_gradients(patch_solution, gradients);
          fe_values[pressure].get_function_values(patch_solution, values);


          for (unsigned int q=0; q<n_q_points; ++q)
            {
              pressure_val += values[q]*values[q]*JxW_values[q];
              for (unsigned int i=0; i<dim; ++i)
                grad_u_val += contract<0,0>(gradients[q][i],gradients[q][i]) * JxW_values[q];
            }
          solu_norm_per_patch += pressure_val + grad_u_val;
        }
    }
  conv_est = std::sqrt(solu_norm_per_patch);
}

// For both p-refinement() and h-refinement() computes the
// error reduction as we named it as patch_convergence_estimator()
template <int dim>
void StokesProblem<dim>::patch_convergence_estimator(const unsigned int cycle,
                                                     SynchronousIterators<std::tuple<DoFHandler_active_cell_iterator,
                                                     std::vector<unsigned int>::iterator>> const &synch_iterator,
                                                     ScratchData &scratch_data, CopyData<dim> &copy_data)
{
  // Silence a warning
  (void) scratch_data;

  DoFHandler_active_cell_iterator cell = std_cxx11::get<0>(synch_iterator.iterators);
  const unsigned int patch_number(*std_cxx11::get<1>(synch_iterator.iterators));
  copy_data.global_cell = cell;
  copy_data.cell_index = patch_number;

  Triangulation<dim> local_triangulation;
  unsigned int level_h_refine(0);
  unsigned int level_p_refine(0);
  std::map<Triangulation_active_cell_iterator, DoFHandler_active_cell_iterator> patch_to_global_tria_map;
  std::vector<DoFHandler_active_cell_iterator> patch = get_patch_around_cell(cell);
  build_triangulation_from_patch (patch, local_triangulation, level_h_refine, level_p_refine, patch_to_global_tria_map);
  hp::DoFHandler<dim> local_dof_handler(local_triangulation);

  set_active_fe_indices (local_dof_handler,patch_to_global_tria_map);
  local_dof_handler.distribute_dofs (fe_collection);

  std::vector<unsigned int> block_component_patch (dim+1, 0);
  block_component_patch[dim]=1;
  DoFRenumbering::component_wise(local_dof_handler, block_component_patch);

  std::vector<types::global_dof_index> dofs_per_block_patch (2);
  DoFTools::count_dofs_per_block(local_dof_handler, dofs_per_block_patch, block_component_patch);

  BlockVector<double> local_solu(dofs_per_block_patch);

  // Here we are trying to project the values of the global vector "solution"
  // into vector "local_solu" which is solution over patch cells corresponding
  // to cell "K".
  DoFHandler_active_cell_iterator act_patch_cell = local_dof_handler.begin_active(),
                                  act_end_patch_cell = local_dof_handler.end();
  for (; act_patch_cell!=act_end_patch_cell; ++act_patch_cell)
    {
      const unsigned int dofs_per_cell = act_patch_cell->get_fe().dofs_per_cell;
      // we check if the corresponding finite element for this cell is not 'FE_Nothing!'
      // and it takes usual finite element.
      if (dofs_per_cell!=0)
        {
          Vector<double> local_solution_values(dofs_per_cell);
          DoFHandler_active_cell_iterator global_cell = patch_to_global_tria_map[act_patch_cell];

          global_cell->get_dof_values(solution,local_solution_values);
          act_patch_cell->set_dof_values(local_solution_values, local_solu);
        }
    }

  BlockVector<double> h_local_solu(local_solu);
  // The local triangulation is not modified by the h-refinement so, we do
  // p-refinement first.
  p_refinement(local_dof_handler, patch_to_global_tria_map, level_p_refine, local_solu);
  patch_solve(local_dof_handler, patch_number, cycle, local_solu, copy_data.p_conv_est_per_cell,
              copy_data.p_workload);

  // Reset the local_dof_handler to what it was before p_refinement was called
  set_active_fe_indices(local_dof_handler,patch_to_global_tria_map);
  local_dof_handler.distribute_dofs(fe_collection);
  DoFRenumbering::component_wise(local_dof_handler, block_component_patch);
  h_refinement(local_triangulation, local_dof_handler, level_h_refine, h_local_solu);
  patch_solve(local_dof_handler, patch_number, cycle, h_local_solu, copy_data.h_conv_est_per_cell,
              copy_data.h_workload);
}

// The marking strategy for cells to be selected for h- or p-refinement
template <int dim>
void StokesProblem<dim>::copy_to_refinement_maps(CopyData<dim> const &copy_data)
{
  h_Conv_Est[copy_data.cell_index] = copy_data.h_conv_est_per_cell;
  p_Conv_Est[copy_data.cell_index] = copy_data.p_conv_est_per_cell;

  const double h_improv = copy_data.h_conv_est_per_cell/est_per_cell[copy_data.cell_index];
  const double p_improv = copy_data.p_conv_est_per_cell/est_per_cell[copy_data.cell_index];

  double h_ratio(h_improv/copy_data.h_workload);
  double p_ratio(p_improv/copy_data.p_workload);

  if (refinement==h_refine)
    p_ratio = 0.;
  if (refinement==p_refine)
    h_ratio = 0.;

  double indicator_per_cell(0.);
  if (h_ratio > p_ratio)
    {
      convergence_est_per_cell[copy_data.cell_index] = h_improv;
      indicator_per_cell = copy_data.h_conv_est_per_cell;
      hp_Conv_Est[copy_data.cell_index] = indicator_per_cell;
      p_ref_map[copy_data.global_cell] = false;

    }
  else
    {
      convergence_est_per_cell[copy_data.cell_index] = p_improv;
      indicator_per_cell = copy_data.p_conv_est_per_cell;
      hp_Conv_Est[copy_data.cell_index] = indicator_per_cell;
      p_ref_map[copy_data.global_cell] = true;
    }

  to_be_sorted.push_back(std::make_pair(indicator_per_cell, copy_data.global_cell));
}

// This function, controls the portion of cells for refinement
// The parameter theta here plays an important role in this step
template <int dim>
void StokesProblem<dim>::mark_cells(const unsigned int cycle, const double theta)
{
  to_be_sorted.clear();
  candidate_cell_set.clear();
  marked_cells.reinit(triangulation.n_active_cells());

  // this vector "convergence_est_per_cell" will be finalized after checking out
  // which h- or p- refinement are going to be chosen for each cell
  convergence_est_per_cell.reinit(triangulation.n_active_cells());

  h_Conv_Est.reinit(triangulation.n_active_cells());
  p_Conv_Est.reinit(triangulation.n_active_cells());
  hp_Conv_Est.reinit(triangulation.n_active_cells());

  // Create the synchronous iterators used by WorkStream
  std::vector<unsigned int> cell_indices(dof_handler.get_tria().n_active_cells(),0);
  for (unsigned int i=0; i<dof_handler.get_tria().n_active_cells(); ++i)
    cell_indices[i] = i;
  DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                  end_cell = dof_handler.end();
  SynchronousIterators<std::tuple<DoFHandler_active_cell_iterator,
                       std::vector<unsigned int>::iterator>> synch_iter(
                         std::tuple<DoFHandler_active_cell_iterator,
                         std::vector<unsigned int>::iterator> (cell, cell_indices.begin()));
  SynchronousIterators<std::tuple<DoFHandler_active_cell_iterator,
                       std::vector<unsigned int>::iterator>> end_synch_iter(
                         std::tuple<DoFHandler_active_cell_iterator,
                         std::vector<unsigned int>::iterator> (end_cell, cell_indices.end()));

  WorkStream::run(synch_iter, end_synch_iter,
                  std_cxx11::bind(&StokesProblem<dim>::patch_convergence_estimator,this, cycle,
                                  std_cxx11::_1, std_cxx11::_2, std_cxx11::_3),
                  std_cxx11::bind(&StokesProblem<dim>::copy_to_refinement_maps,this,
                                  std_cxx11::_1),
                  ScratchData(),
                  CopyData<dim>());

  std::sort(to_be_sorted.begin(), to_be_sorted.end(),
            std_cxx1x::bind(&StokesProblem<dim>::sort_decreasing_order,this,std_cxx1x::_1,std_cxx1x::_2));

  double L2_norm=est_per_cell.l2_norm();
  double sum=0;
  for (unsigned int i=0; i<to_be_sorted.size(); ++i)
    {
      DoFHandler_active_cell_iterator cell_sort = to_be_sorted[i].second;
      sum += (to_be_sorted[i].first)*(to_be_sorted[i].first);

      candidate_cell_set.push_back(cell_sort);

      // if theta is one, refine every cell
      if ((theta<1.0) && (sum >= (theta*(L2_norm))*(theta*(L2_norm))))
        break;
    }

  cell = dof_handler.begin_active();
  for (unsigned int i=0; cell!=end_cell; ++cell,++i)
    {
      typename std::vector<DoFHandler_active_cell_iterator>::iterator  mark_candidate;
      for (mark_candidate=candidate_cell_set.begin();
           mark_candidate!=candidate_cell_set.end(); ++mark_candidate)
        if (cell == *mark_candidate)
          marked_cells[i]=1;
    }
}


template <int dim>
void StokesProblem <dim>::output_results (const unsigned int cycle)
{
  Vector<float> fe_degrees(triangulation.n_active_cells());
  {
    DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                    endc = dof_handler.end();
    for (unsigned int index=0; cell!=endc; ++cell, ++index)
      fe_degrees[index] = fe_collection[cell->active_fe_index()].degree;
  }


  std::vector<std::string> solution_names;
  solution_names.push_back ("x_velocity");
  solution_names.push_back ("y_velocity");
  solution_names.push_back ("pressure");

  std::vector<DataComponentInterpretation::DataComponentInterpretation> data_component_interpretation;
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);
  data_component_interpretation.push_back(DataComponentInterpretation::component_is_scalar);

  DataOut<dim,hp::DoFHandler<dim> > data_out;
  data_out.attach_dof_handler (dof_handler);
  data_out.add_data_vector(solution, solution_names,
                           DataOut<dim,hp::DoFHandler<dim>>::type_dof_data, data_component_interpretation);

  data_out.add_data_vector (marked_cells, "marked_cells");
  data_out.add_data_vector (fe_degrees, "fe_degree");
  data_out.add_data_vector (est_per_cell, "Error_Estimator");
  data_out.add_data_vector (error_per_cell, "Error");
  data_out.add_data_vector (Vect_Pressure_Err, "Pressure_Error");
  data_out.add_data_vector (Vect_grad_Velocity_Err, "Grad_Velocity_Error");

  data_out.add_data_vector (h_Conv_Est, "h_refine_Conv_Est");
  data_out.add_data_vector (p_Conv_Est, "p_refine_Conv_Est");
  data_out.add_data_vector (hp_Conv_Est, "hp_refine_Conv_Est");

  data_out.build_patches ();
  std::string filename = "solution-" +
                         Utilities::int_to_string (cycle, 2) +".vtu";
  std::ofstream output (filename.c_str());
  data_out.write_vtu (output);
}

// This function goes through all cells and tries to do h- or p- refinement
// on all candidade cells which already have been marked for h- or p-refinment

template <int dim>
void StokesProblem<dim>::refine_in_h_p()
{
  bool need_to_h_refine=false;

  typename std::vector<DoFHandler_active_cell_iterator>::iterator  cell_candidate;
  for (cell_candidate=candidate_cell_set.begin(); cell_candidate!=candidate_cell_set.end();
       ++cell_candidate)
    {
      std::vector<DoFHandler_active_cell_iterator> patch_cells =
        get_patch_around_cell (*cell_candidate);

      // Mark the cell for h-refinement
      if (p_ref_map[*cell_candidate]==false)
        {
          need_to_h_refine = true;
          (*cell_candidate)->set_refine_flag();
        }
      // Mark the cell for p-refinement if we haven't reached the maximum
      // polynomial order.
      else if (((*cell_candidate)->active_fe_index()+1) < (fe_collection.size()-1))
        (*cell_candidate)->set_active_fe_index((*cell_candidate)->active_fe_index()+1);
    }

  // Clear the p-refinement map
  p_ref_map.clear();

  if (need_to_h_refine==true)
    triangulation.execute_coarsening_and_refinement();

  bool cell_changed=false;
  do
    {
      cell_changed=false;
      unsigned int count_cell=0;
      DoFHandler_active_cell_iterator cell = dof_handler.begin_active(),
                                      endc = dof_handler.end();
      for (; cell!=endc; ++cell,++count_cell)
        {
          std::vector<DoFHandler_active_cell_iterator> patch_cells = get_patch_around_cell (cell);

          for (unsigned int i=1; i<patch_cells.size(); ++i)
            {

              if (patch_cells[i]->active_fe_index()+1 < (patch_cells[0]->active_fe_index()))
                {
                  patch_cells[i]->set_active_fe_index(patch_cells[0]->active_fe_index()-1);
                  cell_changed=true;
                }
              else if (patch_cells[i]->active_fe_index() > (patch_cells[0]->active_fe_index()+1))
                {
                  patch_cells[0]-> set_active_fe_index (patch_cells[i]->active_fe_index()-1);
                  cell_changed=true;
                }
            }
        }

    }
  while (cell_changed==true);
}


//Explicit initialization
template class StokesProblem<2>;
template class StokesProblem<3>;
