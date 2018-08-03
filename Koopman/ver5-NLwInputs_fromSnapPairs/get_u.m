function u = get_u(t, params)
%get_u: Defines input wrt time (can change if desired)
%   Makes input a sinusoid with specified offset and amplitude

    mean = (params.uub + params.ulb) / 2;
    amplitude = params.uub - mean;
    u = amplitude .* cos(t) + mean;

end