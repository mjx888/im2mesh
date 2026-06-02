function [vert,conn,tria,tnum] = deltri1( node, edge, part )
% deltri1: the simpfied version of deltri2 (revised by Jiexian Ma)
%
% usage:
%	[vert,conn,tria,tnum] = deltri1( node, edge, part );
%
% input:
%   node, edge, part is PSLG data
%
% deltri2 compute a constrained 2-simplex Delaunay triangula-
% tion in the two-dimensional plane.
%   Darren Engwirda : 2017 --
%   Email           : d.engwirda@gmail.com
%   Last updated    : 08/07/2018
%
% copyright (C) by Darren Engwirda


    %---------------------------------------------- extract args
     vert = node;
     conn = edge;
     PSLG = edge;

    %------------------------------------ compute Delaunay tria.
    if (exist('delaunayTriangulation') == +2 )
        dtri = delaunayTriangulation(vert,conn);
        vert = dtri.Points;
        conn = dtri.Constraints;
        tria = dtri.ConnectivityList;
    else
        error('function delaunayTriangulation not exist')
    end

    %------------------------------------ calc. "inside" status!
    tnum = zeros(size(tria,+1),+1) ;

    tmid = vert(tria(:,1),:) ...
         + vert(tria(:,2),:) ...
         + vert(tria(:,3),:) ;
    tmid = tmid / +3.0;

    for ppos = 1 : length(part)
       [stat] = inpoly2( tmid, node, PSLG(part{ppos},:) );
       tnum(stat)  = ppos;
    end
    %------------------------------------ keep "interior" tria's
    tria = tria(tnum>+0,:) ;
    tnum = tnum(tnum>+0,:) ;

    %------------------------------------ flip for correct signs
    area = triarea(vert,tria) ;

    tria(area<0.,:) = tria(area<0.,[1,3,2]) ;

end