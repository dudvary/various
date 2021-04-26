function [] = combineSurfaces(filename1,filename2,outputfile)
% combineSurfaces(file1,file2,outputfile)
% Combine Surfaces
% Input:
% - filename1: path/to/inputSurface1.surf
% - filename2: path/to/inputSurface2.surf
% - outputfile: path/to/combinedSurface.surf

    [Vertices1,Triangles1] = readAmiraSurface(filename1);
    [Vertices2,Triangles2] = readAmiraSurface(filename2);

    Vertices = [Vertices1 Vertices2]; 
    vertexIDOffset = size(Vertices1,2);
    Triangles = [Triangles1 Triangles2+vertexIDOffset]; 
    
    writeAmiraSurface(Vertices,Triangles,outputfile);
    
end