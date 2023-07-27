%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%       POoperation.m                                                     %
%%       PO operations (feed illuminate reflector, reflector to far field, %
%%       etc...                                                            %
%%       Author : Louis Dufour                                             %
%%       Revisions : 0.1 - 16/10/2018                                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

classdef (Abstract) POoperation < matlab.mixin.Copyable
    % POOperation is the Abstract class that all operations must inherit
    % from (such as currents or farfields). Right now, it doesn't do
    % anything at all, but it may be usefull in the future.
    % Authors : Louis Dufour
    % Revision : v0.1.0 - 25/10/2018
end