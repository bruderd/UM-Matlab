function u = gen_fakeInput(t)
%gen_fakeInput: Defines input wrt time (can change if desired)
%   Makes input a sinusoid with specified offset and amplitude
%   Used to generate fake systme data for debugging rest of code

    mean = 0;
    amplitude = 1;
    u = amplitude .* cos(t) + mean;


end

