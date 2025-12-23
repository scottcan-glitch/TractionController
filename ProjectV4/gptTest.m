test_x = x0; test_u = u0;
test_xp = full( xNext(test_x, test_u) );     % integrator function
test_xp_mp = full( xNext_mp(test_x, test_u) );% mpctools wrapper
disp(test_xp'); disp(test_xp_mp');
