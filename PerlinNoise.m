classdef PerlinNoise
    
	properties (SetAccess = private, GetAccess = public)
		
	% Noise dimensionality
		dimension 
			
	% Noise frequencies
		frequencies

	% Intensity of noise frequencies
		amplitudes 		

	% Finite set of gradient directions
		directions 
	
	% Gradients for each noise frequency (selection of <directions>)	
        gradients 		
	
	end
    
    methods
        
		function obj = PerlinNoise( varargin )
			% SYNTAX
			% Obj = PerlinNoise
			% Obj = PerlinNoise( dim )
			% Obj = PerlinNoise( ___, Attribute, Value )

			% DESCRIPTION
			%	Obj = PerlinNoise
			%		Build up a PerlinNoise object with default attributes.
			%
			%	Obj = PerlinNoise( dim )
			%		Build up a PerlinNoise object of dimension <dim>.
			%
			%	Obj = PerlinNoise( ___, Name, Value )
			%		Build up a PerlinNoise object using one or more pairs of 
			%		<Name, Value> noise attributes.
			%
			% ATTRIBUTES ( Name, Value )
			%	( "dir", d ) : 
			%		Specifies the number of spatial directions which will be
			% 		used as gradient vectors. <d> should be a not-null natural 
			%		scalar (i.e. 1,2,3...)
			%
			%	( "freq", vf ) :
			%		Specifies the noise frequency or frequencies. <vf> should 
			%		be a vector (or scalar) of not-null natural numbers. If 
			%		"amp" is specified, both vectors should have the same size.
			%
			%	( "amp", va ) :
			%		Specifies the intensity (amplitude) of noise frequency or 
			%		frequencies. <va> should be a vector (or scalar) of 
			%		not-null positive numbers. If "vf" is specified, both 
			%		vectors should have the same size.
			%
			%	Note: 
			%		- 	Given both vectors <vf> and <va>, the noise frequency 
			%			<vf(i)> will be amplified by  <va(i)>.


		% IMPLEMENTATION

			% Set default values
			obj.dimension 		= 1;
			nDirections 		= 64;
			obj.frequencies 	= 8;
			obj.amplitudes 		= 1;

			% Input processing
			if nargin > 0

				% An odd number of arguments should be given.
				assert( mod(nargin-1,2)==0, ...
					'PerlinNoise(varargin) does not accept an even number of arguments.' );

				% Checking the dimension argument
				assert( isscalar( varargin{1} ), ...
					'PerlinNoise(dim) : <dim> must be an scalar value.' )
				assert( isnumeric( varargin{1} ), ...
					'PerlinNoise(dim) : <dim> must be a number.')
				assert( varargin{1} > 0 , ...
					'PerlinNoise(dim) : <dim> must be greater than 0.')

				if varargin{1} ~= ceil( varargin{1} )

					varargin{1} = ceil( varargin{1} );
					warning( ...
						'PerlinNoise(dim) : <dim> was ceiled up to %i.', ...
						varargin{1} );

				end

				% General attribute name pre-processing
				ordfrom = ones( (nargin - 1)/2 + 1 , 1 );
				ordto = ones( (nargin - 1)/2 + 1 , 1 );
				for i = 2 : 2 : nargin

					% Checking
					assert( isrow( varargin{i} ) && ischar( varargin{i} ), ...
						'PerlinNoise(...,Name(i),...) : <Name(i)> must be ')
					
					% Normalization
					varargin{i} = lower( varargin{i} );

					% Ordering
					tmp = find( strcmp( varargin{i}, ...
						{'dir','freq','amp'} ) ) ;

					assert( ~isempty(tmp), ...
						'PerlinNoise(...,Name(%i),...) : ''%s'' is not a valid attribute name.', i, varargin{i} );
					
					ordfrom( i/2 + 1 ) = i + 1 ;
					ordto( i/2 + 1) = tmp + 1;
					
				end

				% Padding varargin as required for standardization
				tmp = cell( 7 - nargin , 1) ;
				varargin = { varargin{:}, tmp{:} }; 
				
				% varargin reordering
				tmp = cell(1, 5);
				tmp( ordto ) = varargin( ordfrom );
				varargin =tmp;

				% General attribute value checking
				% Number of gradient directions
				if ~isempty( varargin{2} )

					assert( isnumeric( varargin{2} ) && isscalar( varargin{2} ) && varargin{2} > 0, ...
						'PerlinNoise(...,''dir'',d,...) : <d> must be a natural number (scalar) greater than 0.')

				end

				% Noise frequencies and amplitudes relationship
				tmp = cellfun( @(x) numel(x), varargin(3:4) );
				if all(tmp > 0) % Both attributes are specified

					% Check if one, and only one, of the attributes specified 
					% is a scalar:  
					if sum( tmp == 1 ) == 1 

						% Then, convert it in vector of same size with repeated 
						% values ( original varargin{ tmp == 1 }).
						varargin{ tmp == 1 } = repmat( varargin{ tmp == 1 }, 1, tmp( tmp ~= 1 ) );

					else % Check frequency and amplitud vector size.
	
						assert( numel( varargin{3} ) == numel( varargin{4} ),...
							'PerlinNoise(...,''freq'',vf,''amp'',va,...) : <vf> and <va> must be the same size.');

					end

				elseif tmp(2) % Only frequencies are specified 

					% Use default amplitud to fill <va> vector
					varargin{4} = repmat( obj.amplitudes, numel( varargin{3} ), 1 );

				elseif tmp(1) % Only amplitudes are specified

					% Use default frequency to fill <vf> vector
					varargin{3} = repmat( obj.frequencies, numel( varargin{4} ), 1 );

				end

				% Make <vf> and <va> row vectors.
				varargin{3} = reshape( varargin{3}, 1, [] );
				varargin{4} = reshape( varargin{4}, 1, [] );

				% <vf> checking
				assert( isnumeric( varargin{3} ) && all( varargin{3} > 0 ), ...
					'PerlinNoise(..., ''freq'',vf,...) : <vf> must be a positive natural number vector(scalar)') ;

				if ~all( varargin{3} == ceil( varargin{3} ) )

					varargin{3} = ceil( varargin{3} );
					warning('PerlinNoise(..., ''freq'',vf,...) : <vf> was ceiled up to [ %s]', ...
						sprintf('%d ', varargin{3}));

				end

				% <va> checking
				assert( isnumeric( varargin{4} ) && all( varargin{4} > 0 ), ...
					'PerlinNoise(..., ''amp'',va,...) : <va> must be a positive  number vector(scalar)') ;

				% Sorting frequencies - amplitudes : 
				%	Highest frequency should be located at obj.frequencies(1)
				[ varargin{3}, tmp ] = sort( varargin{3}, 'descend' );
				varargin{4} = varargin{4}(tmp);

				% Assignment
				obj.dimension 			= varargin{1};
				if numel( varargin ) > 1 
					if ~isempty( varargin{2} )
						nDirections 		= varargin{2};
					end
					if ~isempty( varargin{3} )
						obj.frequencies 	= varargin{3};
					end
					if ~isempty( varargin{4} )
						obj.amplitudes		= varargin{4};
					end
				end

			end

			% Compute gradient directions
			if obj.dimension > 1
				obj.directions = PerlinNoise.urand( obj.dimension, nDirections );
			else
				obj.directions = rand( nDirections, 1 );
			end
			
			% Compute grid gradients for every frequency
			nF = numel( obj.frequencies );
			obj.gradients = cell( nF, 1 );
			for i = 1 : nF

				obj.gradients{i} = randi(nDirections, obj.frequencies(1).^obj.dimension,1);

			end
            
			
		end
		
		function values = eval(obj, varargin)
			% SYNTAX
			% values = obj.eval( P )
			% values = obj.eval( X0,X1,..XN )
			% values = obj.eval( ___, option, value )

			% DESCRIPTION
			%	values = obj.eval( P )
			%		Evaluates perlin noise for points P. P must be a matrix 
			%		[NxD] where N is the number of points and D its 
			%		dimensionality (which has to be equal to the noise 
			%		dimension).
			%
			%	values = obj.eval( X0,X1,...XN )
			%		Evaluates perlin noise on the N-grid which axis coordinates 
			%		are X0, X1,.. XN. Provided arguments must be equal to noise 
			%		dimension.
			%
			%	values = obj.eval( ___, option, value ) *
			%		Pairs of option-value can be given to modify the evalation 
			%		process. Follow options are available:
			%
			%			- 'normalized' : 
			%				Indicates if noise value should be normalized. 
			%				Default value is True
			%			- 'interp' :
			%				Indicates which interporlation method to use. 
			%				Available methods are: 
			%				'c0' : Zero-continuity, linear interpolation.
			%				'c1' : One-continuity, cubic interpolation.
			%				'c2' : Two-continuity, quintic interporlation.
			%				Default is 'c1'
			
		% IMPLEMENTATION
	
			assert( nargin > 1, 'obj.eval() : No argument was given.');

			% Default parameters
			bNormalized 	= true;
			itpMode			= 'c1';

			% Check <P> and <X0,X1,..XN> arguments
			if obj.dimension == 1

				assert( isvector( varargin{1} ), ...
					'obj.eval(X0) : <X0> must be a vector.');
				
				assert( isnumeric( varargin{1}), ...
					'obj.eval(X0) : <X0> must be numeric.');

				P = reshape( varargin{1}, [], 1) ;
				varargin(1) = [];

			else % For obj.dimension > 1

				assert( ismatrix( varargin{1} ), ...
					'obj.eval( P ) : <P> must be a matrix.') ;

				if size( varargin{1}, 2 ) == obj.dimension

					assert( isnumeric( varargin{1} ), ...
						'obj.eval( P ) : <P> must be numeric.');

					P = varargin{1};
					varargin(1) = [];

				else

					assert( nargin > obj.dimension, ...
						'obj.eval( X0,X1...XN ) : At least %i arguments must be given', obj.dimension );
					
					tmp = cellfun( @(x) isvector(x) && isnumeric(x), varargin( 1 : obj.dimension ) );
					assert( all(tmp), ...
						'obj.eval( X0,X1...XN ) : All <X> must be numeric vector.');

					[ lattice{ 1 : obj.dimension } ] = ndgrid( varargin{ 1 : obj.dimension });
					P = cellfun( @(x) reshape( x, [], 1 ), lattice, 'UniformOutput', false );
					P = cell2mat(P);
					varargin( 1 : obj.dimension ) = [];

				end

			end


			% % Check arguments
			% assert( ( nargin == 2 ) || ( nargin == 1 + obj.dimension ), 'PerlinNoise.eval() method accepts a list of %i-dimensonal points P , or a grid with X1,... X%i vector of axis values',obj.dimension,obj.dimension); 
			% if nargin == 2
			
			% 	assert( isnumeric(varargin{1}) && ismatrix(varargin{1}) && ( size( varargin{1}, 2 ) == obj.dimension ), 'PerlinNoise.eval(P) - P must be a non-empty Nx%i matrix', obj.dimension);
			% 	P = varargin{1};
				
			% elseif nargin == 1 + obj.dimension
				
			% 	for k = 1 : (nargin - 1)
			% 		assert( isnumeric(varargin{k}) && isvector(varargin{k}), 'PerlinNoise.eval( X1,...X%i ) - X%i must be a numeric vector', obj.dimension, k )
			% 	end
				
			% 	[ lattice{ 1: (nargin-1) } ] = ndgrid( varargin{1:(nargin-1)} );
			% 	P = cellfun( @(x) x(:), lattice, 'UniformOutput', false );
			% 	P = cell2mat(P);
				
			% end
			
			% Compute noise for points P
			values = obj.generate( P , @( x ) 1 - x.^2.*(3-2.*x) );
			
			% Normalization
			values = values - min(values(:));
			values = values ./ max(values(:));
			
			% Reshaping
			if ( nargin == 1 + obj.dimension ) && ( obj.dimension > 1 )
				values = reshape( values, cellfun( @(e) numel(e), varargin ) );
			end
			
		end
        
    end
	
	methods ( Access = private )

		function values = generate( obj, P , wf )
			% Computes noise for a given PerlinNoise object and a list of 
			% n-dimensonal points. It gives back a non-normalized vector of 
			% values (scalars) with size equals to the number of points. It 
			% uses wf as interpolation function between corners.

			% Number of points
			nP = size( P , 1 ); 

			% Number of vertices of then n-dimensional cell
			nC = 2 ^ obj.dimension ; 

			% Initialization of output values
			values = zeros ( nP, 1 ); 

			% Relative cell boundary indices
			iB = 0 : ( nC - 1 ) ;
			iB = dec2bin( iB ) ;
			iB = arrayfun( @(x) str2double(x), iB(:) ) ;
			iB = reshape( iB, nC, obj.dimension ) ;

			% Pre-allocate absolute cell boundary indices
			iC = zeros( nP , obj.dimension, nC ) ;

			% Pre-allocate vector indices ( from cell indices to point indices )
			iV = zeros( nP, obj.dimension, nC ) ;

			% Pre-allocate gradient vectors
			iG = zeros( nP, obj.dimension, nC ) ;

			% Loop through noise frequencies
			for f = obj.frequencies 

				% Calculate frequency index
				iF = ( obj.frequencies == f ) ;

				% Point coordinates in the grid
				iP = mod( P, 1 ) .*	f + 1 ; 
				
				% Loop through every cell vertex
				for i = 1 : nC

					% Cell index initialization
					if i == 1
						iC(:,:,i) = floor( iP ) ;
					else
						iC(:,:,i) = mod( iC(:,:,1) + repmat( iB( i, : ), nP, 1 ) - 1, f ) + 1 ;
					end

					% Vector index evaluation
					iV(:,:,i) = iP - iC(:,:,1) - repmat( iB( i, : ), nP, 1 );
					
					% Gradient index evaluation (indices from the grid)
					if obj.dimension > 1
						tmp = mat2cell( iC(:,:,i), nP, ones( 1, obj.dimension ) );
						tmp = sub2ind( repmat( obj.frequencies(1), 1 , obj.dimension ), tmp{:});
					else
						tmp = iC(:,:,i);
					end
					
					% Gradient evaluation
					iG(:,:,i) = obj.directions( obj.gradients{iF}( tmp ), : ) ;

				end

				% Compute dot-product 
				VG = reshape( sum( iV.*iG, 2 ) + 1 , nP, nC ) / 2 ;
				
				% Weight initialization
				W = ones( nP, nC ) ;

				% Loop through dimension to update weights
				for i = obj.dimension : -1 : 1
						
					w = repmat( wf( iV(:,i,1) ), 1, 2.^(obj.dimension-1)) ;
					tmp = iB(:,i) == 0;
					
					W(:,tmp) = W(:,tmp) .* w ;
					W(:,~tmp) = W(:,~tmp) .* ( 1 - w ) ;
					
				end

				% Value update
				values = values + ...
				obj.amplitudes( iF ) .* sum( VG .* W , 2 ) ;

			end

		end

	end
	
    methods ( Static, Access = private)
        
		function u = urand( dim, num )
			% Generates unitary random vectors, uniformly distributed over an 
			% n-sphere.

			% Algorithm description:
			%	urand works sequentially selecting the n-sphere planes. That 
			%	produces some bias produced by the order of selection and 
			%	should be counteracted by a random permutation of the output 
			%	coordinates. 

			switch dim

				case 1
					% Trivial case where unitary random vectors collapse to 
					% a random stream of -1 and 1.
					u = 2.*randi( 2, num, dim ) - 3 ;

				otherwise
					% Initial random coordinates in range -1 to 1.
					u 	= 2.*rand( num, dim ) - 1;

					% Squared remainder  
					r2 	= ones( num, 1 );

					% Loop through coordinate components from highest down to 
					% 3rd.
					for i = dim : -1 : 3

						% Substract squared component value from the remainder
						r2 = r2 - u( :, i ) .^ 2 ;

						% Re-map next component value to ensure final vectors 
						% are normalized.
						u( :, i-1 ) = sqrt(r2) .* u( :, i - 1) ;

					end
					
					% Substract squared component value from the remainder
					r2 = r2 - u( :, 2 ) .^ 2 ;

					% 1st coordinate component is set according the squared 
					% remainder. 
					u( :, 1 ) = sqrt(r2) .* ( 2.* randi( 2, num, 1) - 3 ) ;

				% Shuffling of coordinates to avoid bias produced by sequential 
				% assignment of coordinates.
					iu = repmat( transpose( 1 : num ), 1, dim ) ;
					ju = repmat( dim, num, 1 );
					ju = arrayfun( @(x) randperm(x), ju, 'UniformOutput', false );
					ju = cell2mat( ju );

					u = u( sub2ind( [ num, dim ], iu, ju ) );

			end

		end
		
    end
    
end