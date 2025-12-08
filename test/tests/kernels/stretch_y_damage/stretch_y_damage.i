###############################################################################
#                                stretch_y_GED.i                              #
###############################################################################
# A ``sign of life'' test which makes sure that a simple stretching test      #
# problem runs to completion using the gradient-enhanced damage material      #
###############################################################################

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 1
  ny = 2
  nz = 1
[]

[GlobalParams]
  displacements = "disp_x disp_y disp_z"
  micro_displacement_gradient = "phi_xx phi_xy phi_xz phi_yx phi_yy phi_yz phi_zx phi_zy phi_zz"
  nonlocal_damage = "nonlocal_damage"
  family = LAGRANGE
[]

[Variables]
  [disp_x] []
  [disp_y] []
  [disp_z] []
  [phi_xx] []
  [phi_yy] []
  [phi_zz] []
  [phi_yz] []
  [phi_xz] []
  [phi_xy] []
  [phi_zy] []
  [phi_zx] []
  [phi_yx] []
  [nonlocal_damage] []
[]


[Kernels]
  #Define the internal force balance equations
  [./force_x]
    type = GradientEnhancedDamagedInternalForce
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
    save_in = force_x
  [../]
  [./force_y]
    type = GradientEnhancedDamagedInternalForce
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
    save_in = force_y
  [../]
  [./force_z]
    type = GradientEnhancedDamagedInternalForce
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
    save_in = force_z
  [../]
  #Define the internal couple balance equations
  [./couple_11]
    type = GradientEnhancedDamagedInternalCouple
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
    type = GradientEnhancedDamagedInternalCouple
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
    type = GradientEnhancedDamagedInternalCouple
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
    type = GradientEnhancedDamagedInternalCouple
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
    type = GradientEnhancedDamagedInternalCouple
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
    type = GradientEnhancedDamagedInternalCouple
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
    type = GradientEnhancedDamagedInternalCouple
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
    type = GradientEnhancedDamagedInternalCouple
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
    type = GradientEnhancedDamagedInternalCouple
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
  [GE_Damage]
    type = GradientEnhancedMicromorphicDamage
    variable = nonlocal_damage
    nonlocal_radius = 1.
    use_displaced_mesh = false
  []
[]

[AuxVariables]
  [force_x][]
  [force_y][]
  [force_z][]
  [./test]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./pk2_11]
    order = FIRST
    family = MONOMIAL
  [../]
  [./pk2_22]
    order = FIRST
    family = MONOMIAL
  [../]
  [./pk2_33]
    order = FIRST
    family = MONOMIAL
  [../]
  [./sigma_11]
    order = FIRST
    family = MONOMIAL
  [../]
  [./sigma_22]
    order = FIRST
    family = MONOMIAL
  [../]
  [./sigma_33]
    order = FIRST
    family = MONOMIAL
  [../]
  [./macro_isv]
    order = FIRST
    family = MONOMIAL
  [../]
  [./micro_isv]
    order = FIRST
    family = MONOMIAL
  [../]
  [./micro_gradient_isv_1]
    order = FIRST
    family = MONOMIAL
  [../]
  [./micro_gradient_isv_2]
    order = FIRST
    family = MONOMIAL
  [../]
  [./micro_gradient_isv_3]
    order = FIRST
    family = MONOMIAL
  [../]
  [Fp_22]
    order = FIRST
    family = MONOMIAL
  []
  [omega]
    order = CONSTANT
    family = MONOMIAL
  [] 
[]

[AuxKernels]
  [./test]
    type = MaterialStdVectorAux
    property = PK2
    index = 0
    variable = test
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

[AuxKernels]
  [./macro_isv]
    type = MaterialStdVectorAux
    property = SDVS
    index = 50
    variable = macro_isv
  [../]
[]

[AuxKernels]
  [./micro_isv]
    type = MaterialStdVectorAux
    property = SDVS
    index = 51
    variable = micro_isv
  [../]
[]

[AuxKernels]
  [./micro_gradient_isv_1]
    type = MaterialStdVectorAux
    property = SDVS
    index = 52
    variable = micro_gradient_isv_1
  [../]
[]

[AuxKernels]
  [./micro_gradient_isv_2]
    type = MaterialStdVectorAux
    property = SDVS
    index = 53
    variable = micro_gradient_isv_2
  [../]
[]

[AuxKernels]
  [omega]
    type = MaterialStdVectorAux
    property = ge_damage_statevars
    index = 1 
    variable = omega
  [../]
[]

[AuxKernels]
  [./micro_gradient_isv_3]
    type = MaterialStdVectorAux
    property = SDVS
    index = 54
    variable = micro_gradient_isv_3
  [../]
[]

[BCs]
  active = 'left_x back_z bottom_y top_y'
#  active = 'left_x back_z bottom_y bottom_x top_y top_x'
  [./left_x]
    type = DirichletBC
    #type = PresetBC
    variable = disp_x
    boundary = 'left'
    #boundary = 'left right bottom top front back'
    preset = true
    value = 0
  [../]
  [./back_z]
    type = DirichletBC
    #type = PresetBC
    variable = disp_z
    boundary = 'back'
    #boundary = 'left right bottom top front back'
    preset = true
    value = 0
  [../]
  [./bottom_x]
    type = DirichletBC
    #type = PresetBC
    variable = disp_x
    boundary = 'bottom'
    #boundary = 'left right bottom top front back'
    preset = true
    value = 0
  [../]
  [./bottom_y]
    type = DirichletBC
    #type = PresetBC
    variable = disp_y
    boundary = 'bottom'
    #boundary = 'left right bottom top front back'
    preset = true
    value = 0
  [../]
  [./top_x]
    type     = DirichletBC
    #type     = PresetBC
    variable = disp_x
    boundary = 'top'
    preset = true
    value    = 0
  [../]
  [./top_y]
    #type = DirichletBC
    #type = PresetBC
    type = FunctionDirichletBC
    variable = disp_y
    boundary = 'top'
    #boundary = 'left right bottom top front back'
    preset = true
    function = top_bc
  [../]
[]

[Functions]
  [./top_bc]
    type  = ParsedFunction
    value = 0.2*t
  [../]
[]

[Postprocessors]
  [bot_react_y]
    type = NodalSum
    variable = force_y
    boundary = "top"
  []
  [max_omega]
    type = ElementExtremeValue
    variable = omega
  []
  [max_nl_dmg]
    type = NodalExtremeValue
    variable = nonlocal_damage
  []
  [max_macro_isv]
    type = ElementExtremeValue
    variable = macro_isv
  []
[]

[Materials]
  [linear_elastic_perfect_plastic]
    type = GradientEnhancedDamagedMicromorphicMaterial
    model_name = "LinearElasticityDruckerPragerPlasticity"
    material_fparameters = "2 1.0 2.0
                            2 4.0 5.0
                            2 6.0 7.0
                            2 0. 0.
                            2 0. 0.
                            2 0. 0.
                            2 0. 0.
                            2 0. 0.
                            2 0. 0.
                            2 28. 29.
                            5 31. 32. 33. 34. 35.
                            11 37. 38. 39. 40. 41. 42. 43. 44. 45. 46. 47.
                            2 32. 35.
                            0.5 0.5 0.5 1e-9 1e-9"
    user_material_prop_names = "cu0 Hu cchi0 Hchi cnablachi0 Hnablachi lambda mu eta tau kappa nu sigma tau1 tau2 tau3 tau4 tau5 tau6 tau7 tau8 tau9 tau10 tau11 tauD sigmaD"
    user_material_prop_indices = "1 2 4 5 7 8 28 29 31 32 33 34 35 37 38 39 40 41 42 43 44 45 46 47 49 50"

    number_SDVS = 55
    gradient_enhanced_damage_fparameters = '0.015 0.05 1.0 1'
  []
[]

[UserObjects]
  [reader_element]
    type = PropertyReadFile
    prop_file_name = "stretch_y_heterogeneous_parameters.csv"
    read_type = "element"
    nprop = 24
  []
[]

[Materials]
  [E_nu]
    type = GenericFunctionMaterial
    prop_names = "cu0 Hu cchi0 Hchi cnablachi0 Hnablachi lambda mu eta tau kappa nu sigma tau1 tau2 tau3 tau4 tau5 tau6 tau7 tau8 tau9 tau10 tau11 tauD sigmaD"
    prop_values = "func_cu0 func_Hu func_cchi0 func_Hchi func_cnablachi0 func_Hnablachi func_lambda func_mu func_eta func_tau func_kappa func_nu func_sigma func_tau1 func_tau2 func_tau3 func_tau4 func_tau5 func_tau6 func_tau7 func_tau8 func_tau9 func_tau10 func_tau11 func_tau func_sigma"
    outputs = exodus
    []
[]

[Functions]
  [func_cu0]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "0"
  []
  [func_Hu]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "1"
  []
  [func_cchi0]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "2"
  []
  [func_Hchi]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "3"
  []
  [func_cnablachi0]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "4"
  []
  [func_Hnablachi]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "5"
  []
  [func_lambda]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "6"
  []
  [func_mu]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "7"
  []
  [func_eta]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "8"
  []
  [func_tau]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "9"
  []
  [func_kappa]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "10"
  []
  [func_nu]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "11"
  []
  [func_sigma]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "12"
  []
  [func_tau1]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "13"
  []
  [func_tau2]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "14"
  []
  [func_tau3]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "15"
  []
  [func_tau4]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "16"
  []
  [func_tau5]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "17"
  []
  [func_tau6]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "18"
  []
  [func_tau7]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "19"
  []
  [func_tau8]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "20"
  []
  [func_tau9]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "21"
  []
  [func_tau10]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "22"
  []
  [func_tau11]
    type = PiecewiseConstantFromCSV
    read_prop_user_object = "reader_element"
    read_type = "element"
    column_number = "23"
  []
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
  petsc_options_iname = "-pc_type -pc_factor_mat_solver_package"
  petsc_options_value = "lu       superlu_dist                 "

  line_search = none
  automatic_scaling = true
  l_max_its = 10
  nl_rel_tol = 1e-7
  nl_abs_tol = 1e-7
  nl_max_its = 8 
  start_time = 0.0 
  end_time = 1.0
  dtmin = 1e-12

  dtmax= 0.1
  
  dt = 0.05

  #[TimeStepper]
  #  type = IterationAdaptiveDT
  #  optimal_iterations = 6 
  #  growth_factor=2.0
  #  cutback_factor=0.5
  #  dt = 0.025
  #[] 
[]

[Outputs]
  exodus = true
  perf_graph = true
  print_linear_residuals = true
  #interval = 1
  [log_output]
    type = CSV
  []
#  file_base = 'test'
  [pgraph]
    type = PerfGraphOutput
    execute_on = 'timestep_end final'
    level = 6
  []
[]
