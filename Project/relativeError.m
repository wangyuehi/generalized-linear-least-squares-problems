function [ error ] = relativeError( x, xhat )

    error = norm(x - xhat)/norm(x);


end

