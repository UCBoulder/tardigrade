###############################################################################
#                       test_CSVTimeHistoryDirichletBC.i                      #
###############################################################################
# Test that the CSVTimeHistoryDirichletBC can appropriately apply unique      #
# Dirichlet BCs for unique nodes in a nodeset from multiple CSV files for     #
# multiple degrees of freedom. Time-varying BCs are applied for the "top"     #
# nodeset for u1, u2, u3, and phi_33 DOFs. Time-varying BCs are applied for   #
# all nodes for phi_31 DOFs.                                                  #
###############################################################################

[Mesh]
  type = FileMesh
  displacements = "disp_x disp_y disp_z"
  dim = 3
  file = "mesh.e"
[]

[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
  [./phi_xx]
  [../]
  [./phi_yy]
  [../]
  [./phi_zz]
  [../]
  [./phi_yz]
  [../]
  [./phi_xz]
  [../]
  [./phi_xy]
  [../]
  [./phi_zy]
  [../]
  [./phi_zx]
  [../]
  [./phi_yx]
  [../]
[]

[Kernels]
  #Define the internal force balance equations
  [./force_1]
    type = InternalForce
    component = 0
    dof_num   = 0
    variable  = disp_x

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./force_2]
    type = InternalForce
    component = 1
    dof_num   = 1
    variable  = disp_y

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./force_3]
    type = InternalForce
    component = 2
    dof_num   = 2
    variable  = disp_z

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  #Define the internal couple balance equations
  [./couple_11]
    type = InternalCouple
    component_i = 0
    component_j = 0
    dof_num     = 3
    variable    = phi_xx

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_12]
    type = InternalCouple
    component_i = 0
    component_j = 1
    dof_num     = 4
    variable    = phi_xy

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_13]
    type = InternalCouple
    component_i = 0
    component_j = 2
    dof_num     = 5
    variable    = phi_xz

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_21]
    type = InternalCouple
    component_i = 1
    component_j = 0
    dof_num     = 6
    variable    = phi_yx

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_22]
    type = InternalCouple
    component_i = 1
    component_j = 1
    dof_num     = 7
    variable    = phi_yy

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_23]
    type = InternalCouple
    component_i = 1
    component_j = 2
    dof_num     = 8
    variable    = phi_yz

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_31]
    type = InternalCouple
    component_i = 2
    component_j = 0
    dof_num     = 9
    variable    = phi_zx

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_32]
    type = InternalCouple
    component_i = 2
    component_j = 1
    dof_num     = 10
    variable    = phi_zy

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
  [./couple_33]
    type = InternalCouple
    component_i = 2
    component_j = 2
    dof_num     = 11
    variable    = phi_zz

    #Coupled variables
    u1     = disp_x
    u2     = disp_y
    u3     = disp_z
    phi_11 = phi_xx
    phi_22 = phi_yy
    phi_33 = phi_zz
    phi_23 = phi_yz
    phi_13 = phi_xz
    phi_12 = phi_xy
    phi_32 = phi_zy
    phi_31 = phi_zx
    phi_21 = phi_yx
  [../]
[]


[AuxVariables]
  [./pk2_11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pk2_22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pk2_33]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sigma_11]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sigma_22]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./sigma_33]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./pk2_11]
    type = MaterialStdVectorAux
    property = PK2
    index = 0
    variable = pk2_11
  [../]
[]

[AuxKernels]
  [./pk2_22]
    type = MaterialStdVectorAux
    property = PK2
    index = 4
    variable = pk2_22
  [../]
[]

[AuxKernels]
  [./pk2_33]
    type = MaterialStdVectorAux
    property = PK2
    index = 8
    variable = pk2_33
  [../]
[]

[AuxKernels]
  [./sigma_11]
    type = MaterialStdVectorAux
    property = SIGMA
    index = 0
    variable = sigma_11
  [../]
[]

[AuxKernels]
  [./sigma_22]
    type = MaterialStdVectorAux
    property = SIGMA
    index = 4
    variable = sigma_22
  [../]
[]

[AuxKernels]
  [./sigma_33]
    type = MaterialStdVectorAux
    property = SIGMA
    index = 8
    variable = sigma_33
  [../]
[]

[BCs]
  [./left_x]
    type = DirichletBC
    #type = PresetBC
    variable = disp_x
    boundary = 'left'
    preset = true
    value = 0
  [../]
  [./back_z]
    type = DirichletBC
    #type = PresetBC
    variable = disp_z
    boundary = 'back'
    preset = true
    value = 0
  [../]
  [./bottom_y]
    type = DirichletBC
    #type = PresetBC
    variable = disp_y
    boundary = 'bottom'
    preset = true
    value = 0
  [../]
  [./top_y]
    type = CSVTimeHistoryDirichletBC
    variable = disp_y
    boundary = 'top'
    csv_file = 'top_bc_u2.csv'
  [../]
  [./top_x]
    type = CSVTimeHistoryDirichletBC
    variable = disp_x
    boundary = 'top'
    csv_file = 'top_bc_u1.csv'
  [../]
  [./top_z]
    type = CSVTimeHistoryDirichletBC
    variable = disp_z
    boundary = 'top'
    csv_file = 'top_bc_u3.csv'
  [../]
  [./top_phi_33]
    type = CSVTimeHistoryDirichletBC
    variable = phi_zz
    boundary = 'top'
    csv_file = 'top_bc_phi33.csv'
  [../]
  [./all_phi_31]
    type = CSVTimeHistoryDirichletBC
    variable = phi_zx
    boundary = 'all_nodes'
    csv_file = 'all_nodes_phi31.csv'
  [../]
[]

[Materials]
  [./linear_elastic]
    type = MicromorphicMaterial
    material_fparameters = '2 29.48e3 25.48e3
                            5 1e3 0.4e3 -1.5e3 -1.4e3 -3e3
                            11 0 0 0 0 0 0 10e5 0 0 0 0
                            2 .4e3 -3e3'
    model_name = "LinearElasticity"

    #Coupled variables
    u1     = 'disp_x'
    u2     = 'disp_y'
    u3     = 'disp_z'
    phi_11 = 'phi_xx'
    phi_22 = 'phi_yy'
    phi_33 = 'phi_zz'
    phi_23 = 'phi_yz'
    phi_13 = 'phi_xz'
    phi_12 = 'phi_xy'
    phi_32 = 'phi_zy'
    phi_31 = 'phi_zx'
    phi_21 = 'phi_yx'
  [../]
[]

[Preconditioning]
  [./SMP]
    type = SMP
#    type = FDP
    full = true
  [../]
[]

[Executioner]
  type = Transient

  solve_type = NEWTON
  petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
  petsc_options_value = 'lu       superlu_dist                 '

  line_search = none
  automatic_scaling = true

  nl_rel_tol = 1e-8
  nl_abs_tol = 1e-8
  nl_max_its = 50

  start_time = 0.0
  end_time = 1.0
  dt = 0.1
[]

[Outputs]
  exodus = true
  perf_graph = true
  print_linear_residuals = true
  interval = 1
  [log_output]
    type = CSV
  []
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'timestep_end final'
    level = 2
  []
[]
