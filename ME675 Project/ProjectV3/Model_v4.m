function f = Model_v4()

  [State, StateDeriv, Input, InputDeriv, P] = makeSymbols()
                        [Slip, FrontWheels] = makeKinematics(State, StateDeriv, Input, InputDeriv, P)
                                         Fz = makeNormalForces(StateDeriv, P)
                           WheelFrameForces = makeWheelFrameForces(Slip, P, Fz, FrontWheels);
                 [omega_dot_2, omega_dot_3] = makeDriveTrain(State, Input, P, WheelFrameForces);
                  [v_dot_x, v_dot_y, r_dot] = makeBodyFrameDynamics2(State, StateDeriv,Input,P,WheelFrameForces);

                                          f = [v_dot_x, v_dot_y, r_dot, omega_dot_2, omega_dot_3];
end