
function y = downsample(x, n, phase)

if nargin<2 || nargin>3, print_usage; end

if nargin<3
  phase = 0;
end

if phase > n - 1
  warning('This is incompatible with Matlab (phase = 0:n-1). See octave-forge signal package release notes for details.')
end

if isvector(x)
  y = x(phase + 1:n:end);
else
  y = x(phase + 1:n:end,:);
end

