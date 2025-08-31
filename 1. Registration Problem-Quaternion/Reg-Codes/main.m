%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Registration Algorithm in Robotics using quaternion representation of SO(3),
%   special orthogonal group.
%
%------------------------------------------------------------------------
%   Created by: Nili Abtahi
%
%   Contact: Nili Abtahi (abtahi.niloofar.me@gmail.com)
%
% ========================================================================
%    main function
%
%    Remark. A frame stands for a coordinate system in 3-dimentional Euclidean system
%
% ----------------------------
%    Input:
%
%    *    first_frame_data.txt :  Name of a text file containing position vectors of 3D points in a coordinate system
%
%    *    second_frame_data.txt:  Name of a text file containing position vectors of 3D points in
%                                       another coordinate system
%
%    --------------------------------------------------------------------
%
%    Let A_data be array of  3D position vectors in  first_frame_data.txt and
%        B_data be array of  3D position vectors in  second_frame_data.txt.
%
%    Also, let G be an element of SE(3), special Euclidean, group satisfying B_data = G * A_data .
%    The aim is to find G by least square problem.
%    --------------------------------------------------------------------
%
%    Output:
%        The element G in SE(3) group obtained from least square problem will be printed in terminal/cmd (stdout)
%
%
% ----------------------------
%    Use:
%
%       main  first_frame_data.txt  second_frame_data.txt
%
%   Use for test data:
%       main   ./test/data1-firstFrame.txt     ./test/data1-secondFrame.txt
%
% ==============================================================================




  function   main( varargin )

  clc ;               % clear Matlab command window

  % 0. number of arguments must be two
  if( nargin ~= 2 )
      disp("The input should be as:  main  first_frame_data.txt  second_frame_data.txt")
      error("not enough number of argument")
  end

  first_frame_points  = read_input_file( varargin{1}) ;
  second_frame_points = read_input_file( varargin{2}) ;


  [ SE3_group_element , reg_mean_error ] =  data_registration( first_frame_points , second_frame_points )



  % ========================================================================= %
  function  data_matrix = read_input_file( filename )
            file_id  = fopen( filename , 'r' ) ; % % 1. open calbody.txt file and read its content

            if( file_id == -1 )                  % 1.1. check if file was opened successfully
                error( strcat('read_input_file: cannot open input file:', filename) ) ;
            end

            % % 2. Read file content
            tline = fgetl( file_id ) ;           % % 2.2. get first line of input file: an unimportant string
            tline = fgetl( file_id ) ;           % % 2.2. get second line of input file: number of 3D points
            [N,f,~,nextindex] =  sscanf( tline , '%i ' ) ; % % 2.2.1. extract data from the line

            if( isempty(N) || (N < 1) )          % % 2.3. check if number of points is valid
                error('The second line of input file must be number of 3D points') ;
            end

            data_matrix = zeros(N,3, 'double') ;  % % 2.4. define a Nx3 matrix for storing 3D points
            tline = fgetl( file_id ) ;            % % 2.5. get next line: an unimportant string
            for  i = 1 : N                        % % 2.6. loop for reading 3D points
                 tline = fgetl( file_id ) ;      % % 2.6.1. read a line
                 data_matrix( i , 1:3 ) = sscanf( tline , '%f ,  %f  ,  %f ' ) ;  % % 2.6.2. get point data
            end

            fclose( file_id ) ;                   % % 2.7. close the file
  end % of function read_input_files
  % ========================================================================= %



   end  % main function
