library(ChemmineR)
library(ggplot2)

#todo: support $numbered arg by using paste(unique[i]..., 1:3)...

concatenate_plots = function(sdf_list, filename ='./.html/test.jpg', ...) {
	# Wrapper function makes drawing multiplots in one single image easier.
	#
	# Args:
	#	sdf_list: must be a structure that supports length() method and must have length > 1
	# 	filename: filename of multiplot. Include format as .pdf, .jpg, etc in the filename arg
	# 	...: further arguments will be passed to instances of draw_sdf. See ?draw_sdf for more info
	#
	# Returns:
	#	an instance of a ggplot object created with multiplot which contains multiple ggplots
	# Requires:
	#	grid, gridExtra

	require(grid)
	require(gridExtra)
	
	image_list = list();

	for (i in 1:length(sdf_list)) {
		image_list[[i]] = draw_sdf(sdf_list[[i]], filename=NULL, ...)
		} # end iteration through image list

	if(!is.null(filename)) {
		ggsave(filename, do.call(arrangeGrob, image_list))
		}
	} # end function concatenate_plots

get_default_size_dict = function() {
	return(c( 	'default'=1,
			'C'=0,
			'H'=0));
	}
	# todo: where is this used? for what?

get_default_edge_color_dict = function() {
	return(c( 	'default'='black',
			'H'='delete'));
	} # todo deprecated... will be reused in future

get_default_text_color_dict = function() {
	return(c(	'S'='green',
			'F'='gold',
			'N'='blue',
			'O'='red',
			'C'='black',
			'H'='white',
			'default'='black'));
	} # end function to return default color dicts for text coloring

convert_fancy_edges = function(sdf, atom_names, edge_color_dict) {
	# Converts sdf data about bond strength and position into a matrix of single edges to be drawn.
	#
	# Args:
	#	sdf: an sdf to draw
	#	atom_names: list of strings of atoms within the sdf
	#	edge_color_dict: deprecated. a mapping which describes what color each edge should be. Will be replaced
	#		when mcs support is added
	# Returns:
	#	a list of colors which will be columnwise bound to a list of edges
	# Requires:
	#	ChemmineR
	priority_count = 4; # col1 = index of atom 1, col2 = index of atom 2, col3 = default value, col4 = if bond between atoms is defined as having special coloration
	color_text_mat = cbind(atom_names[bondblock(sdf)[,1]], atom_names[bondblock(sdf)[,2]]);

	priority_mat = matrix(data=0, nrow=dim(bondblock(sdf))[1], ncol=priority_count);
	
	priority_mat[,1] = match(color_text_mat[,1], names(edge_color_dict));
	priority_mat[,2] = match(color_text_mat[,2], names(edge_color_dict));
	priority_mat[,3] = match('default', names(edge_color_dict));
	priority_mat[,4] = match(paste(color_text_mat[,1], color_text_mat[,2], sep='-'), names(edge_color_dict));
	priority_mat[is.na(priority_mat)] = 0; # replace all NAs with default value
	priority_mat = apply(priority_mat, 1, max)

	return(edge_color_dict[priority_mat])
	}

convert_nodes_to_text_with_color = function(node_frame, atom_names, numbered, node_size_dict, node_color_dict) {
	# Takes a list of nodes, determines where background whitespace needs to be drawn, determines where to draw text.
	#
	# Args:
	# 	node_frame: a data frame which contains x y data for each atom
	#	atom_names: list of strings describing atom names
	# 	numbered: boolean; if true, the ith unique atom with the same name will have the number i appended to its text
	#	node_size_dict: deprecated, will be removed
	#	node_color_dict: deprecated, will be removed
	frame_builder = node_frame;
	size_col = node_size_dict[atom_names];
	size_col[is.na(size_col)] = node_size_dict['default'];
	frame_builder = cbind(frame_builder, size_col)

	toprint = atom_names;
	if (numbered) {
		unique_list = unique(toprint)
		for (i in 1:length(unique_list)) {
			result = toprint[toprint==unique_list[i]]
			result = paste(result, 1:length(result), sep='-')
			toprint[toprint==unique_list[i]] = result
			} # end iteration through list of unique characters using i
		}
	toprint[size_col==0] = '';
	frame_builder = cbind(frame_builder, toprint);

	color_col = node_color_dict[atom_names];
	color_col[is.na(color_col)] = node_color_dict['default'];

	frame_builder = cbind(frame_builder, color_col);
	
	horiz_text_spacer = 0.001; # constant value for now, todo: make into assignable variablea
	vert_text_spacer = 0.1;

	# create space for xmin, xmax, ymin, ymax for geom_rect to clear space
	frame_builder = cbind(	frame_builder, 
				frame_builder['C1'] - horiz_text_spacer*nchar(frame_builder['toprint']), # xmin
				frame_builder['C1'] + horiz_text_spacer*nchar(frame_builder['toprint']), # xmax
				frame_builder['C2'] - vert_text_spacer, # ymin, does not scale with length of text to print
				frame_builder['C2'] + vert_text_spacer) # ymax, does not scale with length of text to print

	colnames(frame_builder) = c('C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', 'C8', 'C9');
	frame_builder = frame_builder[frame_builder[,3]!=0, ];

	return(data.frame(frame_builder));	
	} # end function to generate a matrix containing X, Y, text, and color data

convert_strong_edge_colors = function(bond_str_vec, incolors) {
	# Converts list of strengths and colors to repeated list of colors of bonds
	#
	# Args:
	#	bond_str_vec: vector of values > 1, describing the strength of the ith bond
	# 	incolors: list of strings of colors to draw, describing the color of the ith edge
	#
	# Returns:
	#	list of colors, similar to incolors but with elements repeated
	#		as many times as corresponding to bond strength

	# could just write cbind(color_vec[j] for i in 1:bond str,for j in color vec)
	if(length(bond_str_vec) != length(incolors)) {
		print(length(bond_str_vec));
		print(length(incolors));
		stop('Error: bond_str_vec and incolors have differing lengths.');
		} # end safety clause to ensure that incolors and bond_str_vec have same length

	outcolor_vec = rep(0,sum(bond_str_vec));
	
	i = 1; j = 1;
	while (i < length(outcolor_vec)) {
		for (k in 1:bond_str_vec[j]) {
			# use i to track which index of outcolor_vec we're on
			# use j to track which index of incolors/bond_str_vec we're on
			# use k to track which bond in bond_str_vec we're on
			outcolor_vec[i] = incolors[j];
			i = i +1;
			}
		j = j + 1;
		#i = i + 1;
		} # end iteration through incolor vec using i

	return(outcolor_vec)	
	} # end function convert_strong_edge_colors

convert_strong_bonds = function(strong_bond_mat, bond_dist_offset = 0.25) {
	# Converts bonds with strength exceeding 1 into a new pair of bonds
	#
	# Args:
	#	strong_bond_mat: matrix describing the locations where a bond originates from
	# 		and terminates as [x1,y1,x2,y2] and bond strength
	#	bond_dist_offset: describes how far apart to draw two+ bonds. 
	#		units are ggplot units, not ChemmineR units
	# Returns:
	#	split_bond_mat: a matrix which contains agnostic information about
	# 		origins and termini of multiple edges, regardless of how strong the 
	#		original bond was

	# todo: change to rbind, sacrifice speed for clarity
	if ((dim(strong_bond_mat))[1] <= 0) {
		stop('Error: No bonds found in convert_strong_bonds.');
		return(NULL); # can rbind NULL into a matrix for no effect
		# todo: expand to allow for multiple images to be catted together
		} # end check to ensure matrix dimension is correct
	
	color_text = strong_bond_mat[,6];
	strong_bond_mat = cbind(strong_bond_mat[,1:5], 0);
	class(strong_bond_mat) = 'numeric';
	split_bond_mat = matrix(data=0, nrow=sum(strong_bond_mat[,5]), ncol=dim(strong_bond_mat)[2]);

	i = 1; j = 1; # i iterates through rows of original strong_bond_mat that describes original multibond edges
	# j iterates through rows of new split_bond_mat that describes each individual edge's start + end in x,y coordinates
	offset_vec = c(1,-1,1,-1);
	while (i <= dim(strong_bond_mat)[1]) {
		dxdy = c(strong_bond_mat[i,1] - strong_bond_mat[i,3], strong_bond_mat[i,2] - strong_bond_mat[i,4]);
		scaled_dxdy = rev((dxdy/sqrt(sum(dxdy^2)))) * (bond_dist_offset);
		dp1 = strong_bond_mat[i,1:4] + scaled_dxdy *  1*offset_vec;
		dp2 = strong_bond_mat[i,1:4] + scaled_dxdy * -1*offset_vec;
		
		split_bond_mat[j,1:4] = dp1;
		split_bond_mat[j,5:6] = c(1,color_text[i]);
		j = j+1;

		split_bond_mat[j,1:4] = dp2;
		split_bond_mat[j,5:6] = c(1,color_text[i]);
		j = j+1;

		if (strong_bond_mat[i,5] > 2) {
			split_bond_mat[j,1:4] = strong_bond_mat[i,1:4]; j = j + 1;
			# reinsert old edge back into matrix
			} # end check: if bond strength exceeds 2
		i = i + 1; # i should be used to iterate through rows of strong_bond_mat
		} # end iteration through rows of strong_bond_mat
	return(split_bond_mat);
	} # end function convert_strong_bonds

convert_fmcs_results = function(...) {
	# todo need to ask Kevin how to provide subfunction arguments
	# from main body function without having 20+ args
	require(fmcsR)
	
	result = fmcs(...)
	# for now, just return fmcs results

	return(result);
	} # end function convert_fmcs_results

draw_sdf = function(	sdf,
			filename = './sdf_draw_test.jpg', # todo: fix this before submission to repository
			node_size_dict = get_default_size_dict(),		# todo: deprecated, remove
			node_color_dict = get_default_text_color_dict(),	# list mapping that pairs key elements with value colors
			edge_color_dict = get_default_edge_color_dict(),	# todo: deprecated, will be recruited in the future to draw e.g. red edges for MCS, for now just tolerate poor implementation
			numbered= FALSE,	# todo: not yet supported well... change from atomnames to paste(1...count + text)
			bgcolor = rgb(1,1,1,1),	#full red, full green, full blue, fully opaque
			bond_dist_offset = 0.05,# how many graphical units to offset nonsingular bonds
			alpha_text = 1.0,	# transparency of text
			alpha_edge = 1.0,	# transparency of edges/bonds
			font_size = 7, 		# note: special behavior when using multiplot... font sizes ARE RESCALED for better readability when multiplotting... todo: enable special behavior
			bg_raster = NULL, 	# user can specify background imaage to bind. When multiploting, imposed_raster must still be a single raster. use readPNG to prepare a raster
			erase_node_background=TRUE
			) {
	# Draws an sdf
	#
	# Args:
	# 	sdf: can be a list of sdfs, in which case concatenate_plots is recruited.
	#		if a single sdf is passed, will attempt to draw that sdf using ggplot.
	#	filename: filename or filepath + filename to draw to. Include format in filename
	#		e.g. "filename = './draw_sdf_test/mymo.pdf'". Supports .jpg and .png
	#	node_size_dict: deprecated
	# 	edge_color_dict: deprecated, will be used when mcs support is added
	#	numbered: boolean describing whether atoms should have numeric values counting their procession
	#	bgcolor: background color of ggplot. Can be any valid color, like "'red'" 
	#		or "rgb(1,1,1)" or "rgb(0.5,0,0.25,1)"
	# 	bond_dist_offset: describes how far apart, in ggplot units, double+ bonds should be spaced
	#	alpha_text: the alpha/transparency of text to be drawn
	#	alpha_edge: the alpha/transparenc of edges to be drawn
	#	font_size: text font size. Note: interacts poorly with multiplots
	#	bg_raster: a raster that has been created using readPNG. Will be set as background
	#	erase_node_backgruond: boolean. If true, draws *whitespace* in a square behind each text box
	#
	# Returns:
	#	a ggplot object
	#	will also save a file at $filename if filename is given
	#	note that if a print() command is even implicitly called on a ggplot object, R will attempt
	#	to draw to your default graphics device
	#
	# Dependencies:
	#	ChemmineR, grid, gridExtra
	#
	# Example of usage:
	# 	library(ChemmineR)
	#	data(sdfsample)
	# 	source('draw_sdf.r'); draw_sdf(sdfsample[[30]])
	
	if (length(sdf) > 1) {
		return(concatenate_plots(sdf[1:length(sdf)], filename, disable_drawing = TRUE, node_size_dict, node_color_dict,
					edge_color_dict, numbered, bgcolor, bond_dist_offset, alpha_text, alpha_edge,
					font_size, bg_raster, erase_node_background))
		}

	# step 1: prepare frames describing the data we're drawing

	node_frame = data.frame(atomblock(sdf))[,1:2];
	colnames(node_frame) = c('C1', 'C2');
	edge_frame = data.frame(bondblock(sdf))[,1:3];
	
	atom_names = gsub('[^[:alpha:]]', '', rownames(atomblock(sdf))); # delete non-alphabetical characters from char array of atom names	
	labeled_node_frame = cbind(node_frame, atom_names); # , rownames(node_frame));
	if (numbered) labeled_node_frame[,3] = rownames(node_frame);
	colnames(labeled_node_frame) = c('C1', 'C2', 'atom_names');
	edge_colors = convert_fancy_edges(sdf, atom_names, edge_color_dict);

	segment_frame = cbind(	node_frame[edge_frame[,1],1],
				node_frame[edge_frame[,1],2],
				node_frame[edge_frame[,2],1],
				node_frame[edge_frame[,2],2],
				edge_frame[,3],
				0)
	colnames(segment_frame) = c('X1' ,'X2', 'X3', 'X4', 'X5', 'X6');
	rownames(segment_frame) = NULL;

	# step 1.5: do convoluted reshapings of matrices involving double+ bonds or different colors
	# because for some of these functions, including both text and numeric data in the same matrix = NOT OK
	strong_bonds = segment_frame[segment_frame[,5]>1,]; 
	strong_edge_colors = edge_colors[segment_frame[,5]>1];
	strong_edge_weights = segment_frame[segment_frame[,5]>1,5];
	
	edge_colors = edge_colors[!segment_frame[,5]>1]; 
	edge_colors = c(edge_colors,convert_strong_edge_colors(strong_edge_weights, strong_edge_colors));
	
	segment_frame = segment_frame[!segment_frame[,5]>1,]; 
	if (length(strong_edge_colors) > 0) {
		segment_frame = rbind(segment_frame, convert_strong_bonds(strong_bonds, bond_dist_offset));
		}
	segment_frame = segment_frame[edge_colors!= 'delete', ]
	edge_colors = edge_colors[edge_colors!= 'delete'];
	segment_frame = data.frame(segment_frame);

	# step 2: begin drawing
	plot_target = ggplot(data=data.frame(0), aes());

	# if user specifies a raster to append
	if (!is.null(bg_raster)) {
		require(png)
		require(grid)
		g = rasterGrob(bg_raster, interpolate = TRUE)
		plot_target = plot_target + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
		}

	plot_target = plot_target + geom_segment(alpha = alpha_edge, data=segment_frame, aes(x=X1,xend=X3, y=X2, yend=X4), color=edge_colors);

	recttext_frame = convert_nodes_to_text_with_color(node_frame, atom_names, numbered, node_size_dict, node_color_dict);
	if (erase_node_background) {
		plot_target = plot_target + geom_rect(data=recttext_frame, aes(xmin=C6, xmax=C7, ymin=C8, ymax=C9), color=bgcolor, fill=bgcolor);
		}
	plot_target = plot_target + geom_text(alpha=alpha_text, data=recttext_frame, aes(x=C1, y=C2, size=C3, label=C4), color = recttext_frame[,'C5'], size=font_size);
	# completed drawing of all doodads

	# step 3: prep for and then export as image

	# make ui beautiful
	#plot_target = plot_target + theme_bw();
	plot_target = plot_target + theme(panel.background= element_rect(fill=bgcolor, color='black'));
	plot_target = plot_target + theme(legend.position='none')
	plot_target = plot_target + theme(axis.title.x = element_blank())
	plot_target = plot_target + theme(axis.title.y = element_blank())
	
	if (!is.null(filename)) {
		ggsave(filename=filename, plot=plot_target);
		}

	return(plot_target);
	} # end function draw_sdf
