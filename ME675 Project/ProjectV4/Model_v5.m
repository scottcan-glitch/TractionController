function Model = Model_v5(car)
disp('starting:')
     [State, StateDeriv, Input, InputDeriv] = makeSymbols2()
        Model.symbols.State = State;
        Model.symbols.StateDeriv = StateDeriv;
        Model.symbols.Input = Input;
        Model.symbols.InputDeriv = InputDeriv;
                    disp('makeKinematics:')
                    [Slip, FrontWheels, wc] = makeKinematics(State, StateDeriv, Input, InputDeriv, car)
                                         disp('makeNormalForces:')
                                         Fz = makeNormalForces(StateDeriv, car)
                                         disp('makeWheelFrameForces:')
                                         WheelFrameForces = makeWheelFrameForces2(Slip, car, Fz, FrontWheels)
                                         disp('makeDriveTrain:')
                                         [omega_dot_2, omega_dot_3] = makeDriveTrain2(State, Input, car, WheelFrameForces)
                                         disp('makeBodyFrameDynamics2:')
                                         [v_dot_x, v_dot_y, r_dot] = makeBodyFrameDynamics3(car, State, StateDeriv, Input, WheelFrameForces, wc)
                                         symvar(omega_dot_2)
                                         symvar(omega_dot_3)
                                         symvar(v_dot_x)
                                         symvar(v_dot_y)
                                         symvar(r_dot)

                                         disp('Next, f:')
                                         % StateDeriv.vector
                                         % size(StateDeriv.vector)
                                         f = StateDeriv.vector' - [v_dot_x; v_dot_y; r_dot; omega_dot_2; omega_dot_3]

                                         f2 = subs(f,'delta_dot',0);    %Very slow steer rate assumption is valid
                                         %f = DAE2ODE(StateDeriv, v_dot_x, v_dot_y, r_dot, omega_dot_2, omega_dot_3)
                                          % Return struct
                                          Model.f = f2
                                          %Model.A = getJacobian(f,State)
                                         % Model.B = getJacobian(f,Input)

end