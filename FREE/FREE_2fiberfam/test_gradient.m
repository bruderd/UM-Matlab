function [ grad, num_grad ] = test_gradient( f, g, x, Prob, use_qe )

test = f( x, Prob );
len_f = length( test );
len_x = length( x );
h = 1e-6;
if( nargin < 5 )
  use_qe = false;
end
qe = [];
if( use_qe )
  qe = rand( len_f, 1 ) <= 0.7;
end

if( use_qe )
  num_grad = zeros( sum( qe ), len_x );
else
  num_grad = zeros( len_f, len_x );
end
for k = 1:len_x,
  hvec = zeros( size(x) );
  hvec(k) = h;
  fph = f( x + hvec, Prob );
  fmh = f( x - hvec, Prob );
  if( any( isnan( fph ) ) || any( isnan( fmh ) ) )
    error( [ 'NaN found at k = ' num2str( k ) ] );
  end
  if( use_qe )
    num_grad(:,k) = ( fph(qe) - fmh(qe) ) / ( 2 * h );
  else
    num_grad(:,k) = ( fph - fmh ) / ( 2 * h );
  end    
end

if( use_qe )
  grad = g( x, Prob, qe );
else
  grad = g( x, Prob );
end

% dan's hack 11/28
num_grad = num_grad';

fprintf( 1, 'Maximum error = %d\n', max( abs( grad(:) - num_grad(:) ) ) );
% if( max( abs( grad(:) - num_grad(:) ) ) > 1e-5 )
%   imshow( abs( grad - num_grad ) > 1e-5 );
% end
