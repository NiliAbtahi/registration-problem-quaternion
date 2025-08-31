% ==============================================================================
%    Registration Algorithm in Robotics using quaternion representation of SO(3),
%    special orthogonal group.
%
%    This function takes two set of 3D points \{ a_ i \} and \{ A_i \}
%    and finds an Euclidean transformation F [ R, t ] \in SE(3)
%    that solve the least square problem:  A_i ~ R a_i + t.
%    The algorithm and equations are discussed in pdf document.
%
%    This function is not called by user directly.  Instead, the main function calls it.
%
%    Input:
%              domain_data  :  the set of 3D point { a_ i } in the least square problem A_i ~ R a_i + t.
%              range_data   :  the set of 3D point { A_ i } in the least square problem A_i ~ R a_i + t.
%
%   Output:
%             SE3_group_element : the matrix F mentioned above.
%             reg_mean_error    : mean Euclidean error of the registration.
%
%    Author:   Nili Abtahi
%    Email  :   abtahi.niloofar.me@gmail.com
%
%
% ==============================================================================

  function  [ SE3_group_element , reg_mean_error ] =  data_registration( domain_data , range_data )

           % % 0. check data validity
           if( size(domain_data , 1) ~= size(range_data , 1))
               error("Number of points in two coordinate systems should be the same")
           end

           if( (size(domain_data , 2) ~= 3) || (size(range_data , 2) ~= 3) )
               error("The points should be in 3D space, i.e. size_of_point_array == Nx3")
           end

          % % 1. get mean of domain and range data
          domain_data_mean = mean( domain_data )  ;
          range_data_mean = mean( range_data ) ;

          % % 2. define deviation vectors, i.e. data - mean( data )
          domain_dev = domain_data - domain_data_mean ;
          range_dev   = range_data - range_data_mean ;

          % % 3. form matrix M = sum_{i} P^T_i Q_i
          num_points = size( domain_data , 1 ) ;
          M = zeros( 4 , 4 , 'double' ) ;

          for  i = 1 : num_points

                 % form P matrix
                 P = [ 0.0  , -domain_dev(i, 1:3)  ;
                          domain_dev(i, 1) , 0.0 , domain_dev(i, 3) , -domain_dev(i, 2)  ;
                          domain_dev(i, 2) , -domain_dev(i, 3) , 0.0 , domain_dev(i, 1) ;
                          domain_dev(i, 3) , domain_dev(i, 2) , -domain_dev(i, 1) , 0.0 ] ;

                 % form Q matrix
                 Q = [ 0.0 , -range_dev( i , 1:3) ;
                           range_dev( i , 1) , 0.0 , -range_dev( i , 3) ,  range_dev( i , 2 ) ;
                           range_dev( i , 2) , range_dev( i , 3)  , 0.0  , -range_dev( i , 1)  ;
                           range_dev( i , 3) , -range_dev( i , 2) , range_dev( i , 1) , 0.0 ] ;

                 % % add   transpose(P) * Q to M
                 M = M + transpose(P) * Q ;
          end

          % % 4. get eigenvalue and eigenvector of M
          [ vec , val ] = eig( M ) ;

          % % 5. the unit quaternion associated to rotation matrix is
          % %     the eigenvector with maximum eigenvalue.

          % % 5.1. find maximum eigenvalue
          max_id = 1 ;
          max_value = val(1,1) ;
          for i = 2:4
                if( val(i,i) > max_value )
                     max_id = i ;
                     max_value = val(i,i) ;
                end
          end


          % % 5.2. choose eigenvector with max_eigen_value as desired quaternion.
          q = vec(:, max_id )  ;

           % % 6. form rotation matrix
           R = zeros(3,3) ;

           R(1,1) = 1.0 - 2.0 * (   q(3)^2 +  q(4)^2 ) ;
           R(2,2) = 1.0 - 2.0 * (   q(2)^2 +  q(4)^2 ) ;
           R(3,3) = 1.0 - 2.0 * (   q(2)^2 +  q(3)^2 ) ;

           R(1,2) = 2.0 * ( q(2) * q(3)  - q(1) * q(4)   ) ;
           R(1,3) = 2.0 * ( q(2) * q(4) + q(1) * q(3)  ) ;
           R(2,1) = 2.0 * ( q(2) * q(3) + q(1) * q(4)  ) ;
           R(2,3) = 2.0 * ( q(3) * q(4)  -  q(1) * q(2)  ) ;
           R(3,1) = 2.0 * ( q(2) * q(4)  -  q(1) * q(3)  ) ;
           R(3,2) = 2.0 * ( q(3) * q(4) +  q(1) * q(2) ) ;


           % % 7. get translation subgroup of SE3_element
           t = range_data_mean' - R * domain_data_mean' ;

          % % 8. form SE3_element
          SE3_group_element = eye(4) ;
          SE3_group_element(1:3 , 1:3) = R ;
          SE3_group_element(1:3 , 4   ) = t ;

          % % 9. get mean error of data_registration
          % % mean Euclidean norm of error = ( 1 / N ) sqrt(  sum_{i=1}^N  norm( error_i )^2 )

          reg_mean_error = 0.0 ;

          for  i = 1 : num_points
                 v_domain = domain_data(i , 1:3 ) ;
                 v_range   =  range_data(  i , 1:3 ) ;

                 err_vec =  R * v_domain' + t  - v_range' ;
                 reg_mean_error = err_vec' * err_vec ;
          end

         reg_mean_error = sqrt( reg_mean_error ) / num_points ;

         disp( sprintf("\n\n The SE(3) group element relating the two coordinate systems is in projective representation:"))
         disp( SE3_group_element )
         disp( sprintf("\n\n Mean error of registration algorithm = %15.12f \n" , reg_mean_error))

  end
