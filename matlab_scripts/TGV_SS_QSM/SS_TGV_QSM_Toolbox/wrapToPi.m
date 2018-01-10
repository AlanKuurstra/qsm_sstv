function [ res ] = wrapToPi( phs )
%WRAPTOPI Summary of this function goes here
%   Detailed explanation goes here

res = angle(exp( 1i .* phs ));

end

