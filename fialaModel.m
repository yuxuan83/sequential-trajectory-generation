function F_y = fialaModel(alpha, C, F_z, mu)
    idx = (abs(alpha) < atan(3*mu*F_z/C));
    F_y = zeros(size(alpha));

    F_y(idx) = (-C*tan(alpha(idx)) + C^2/(3*mu*F_z)*abs(tan(alpha(idx))).*tan(alpha(idx)) - C^3/(27*mu^2*F_z^2)*tan(alpha(idx)).^3);
    F_y(~idx) = -mu*F_z*sign(alpha(~idx));
end