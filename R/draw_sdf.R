
library(ggplot2)
library(grid)
library(gridExtra)

default_node_policy = function() {
	return(c(	'N'='blue',
			'O'='red',
			'F'='gold',
			'S'='green',
			'default'='black',
			'H'='skip',
			'fmcs'='magenta',
			'C'='skip'
			))
	}

default_edge_policy = function() {
	return(c(	'default'='black',
			'H'='skip',
			'fmcs'='magenta'
			))
	
	}

concatenate_plots = function(sdf_list, filename=NULL, ...) {

	image_list = list()
	for (i in 1:length(sdf_list)) {
		image_list[[i]] = draw_sdf(sdf_list[[i]], filename=NULL,...)
		}

	image_list = do.call(arrangeGrob, image_list)
	if (!is.null(filename)) {
		#ggsave(filename, do.call(arrangeGrob, image_list))
		ggsave(filename, image_list)
		}
	#image_list = do.call(arrangeGrob, image_list)
	return(image_list)
	# todo figure out how to return concatenated images
	} # end function concatenate_plots

handle_raster = function(plot_target, raster) {
	require(png)
	require(grid)

	# check if it's a filename
	if (typeof(raster) == 'character') {
		raster = readPNG(raster)
		} 

	g = rasterGrob(raster, interpolate=TRUE)
	plot_target = plot_target + annotation_custom(g, xmin=-Inf, xmax=Inf, ymin=-Inf, ymax=Inf)
	return(plot_target)
	}

handle_text = function(sdf, plot_target, alpha_node = 1.0, numbered = FALSE, font_size = 7, node_vertical_offset = 0, node_background_color = FALSE, node_policy = default_node_policy(), fmcs_result = NULL) {
	# handles printing of text to a plot_target
	# Args:
	#	sdf: from main draw_sdf function
	#	plot_target: must be a ggplot2 object
	#	alpha_node: defines transparency of nodes
	#	numbered: default FALSE, causes all text to be printed as
	#		without numbering. if TRUE or 1, numbers will be
	#		appended which count nth atom of type i.
	#		if TRUE but not 1, numbers will be appended which
	#		count nth atom
	#	font_size: defines absolute font size of text
	#	node_vertical_offset: permits user to shift text upward off of atom points in space
	#	node_policy: color and drawing policy for nodes. See default_node_policy()
	#	fmcs_result: if non-null, will attempt to highlight text which corresponds to 
	#		MCS in sdf. If non-null, node_policy MUST support node_policy['fmcs']	

	# handle atom name preparation
	node_names = gsub('[^[:alpha:]]', '',rownames(atomblock(sdf)) )
	toprint = node_names
	if (numbered) {
		# default behavior if numbered is evaluatable to TRUE
		unique_text = unique(node_names)
		for (i in 1:length(unique_text)) {
			atom_count = sum(node_names==unique_text[i])
			toprint[node_names==unique_text[i]] = paste(node_names[node_names==unique_text[i]], 1:atom_count, sep='')
			} # end iteration through unique text names using i
		if (numbered==2) {
			toprint = rownames(atomblock(sdf))
			}
		} # end if: numbered is not FALSE
	fmcs_labels = toprint[fmcs_result]

	# handle sieving out of skippable atoms
	node_frame = atomblock(sdf)[,1:2]
	color_list = node_policy[node_names]
	color_list[is.na(color_list)] = node_policy['default']
	# delete skippable atoms from toprint, node_names, node_frame, color_list
	toprint 	= toprint[	color_list != 'skip']
	node_names 	= node_names[	color_list != 'skip']
	node_frame	= node_frame[	color_list != 'skip', ]
	color_list	= color_list[	color_list != 'skip']

	# draw atoms as text
	node_frame[,2] = node_frame[,2] + node_vertical_offset
	node_frame = data.frame(node_frame)
	
	# add coloration editing here
	if (node_background_color != FALSE) {
		# presumably it's a color...
		plot_target = plot_target + geom_point(data = node_frame, aes(x=C1, y = C2), color = node_background_color, size = font_size + 1)
		}

	plot_target = plot_target + geom_text(alpha = alpha_node, data=node_frame, aes(x=C1, y=C2), label = toprint, color=color_list, size = font_size)

	if (!is.null(fmcs_result)) {
		fmcs_frame = atomblock(sdf)[fmcs_result,1:2]
		fmcs_frame[,2] = fmcs_frame[,2] + node_vertical_offset
		fmcs_frame = data.frame(fmcs_frame)
		plot_target = plot_target + geom_text(alpha=alpha_node, data=fmcs_frame, aes(x=C1, y=C2), label=fmcs_labels, color=node_policy['fmcs'], size = font_size)

		} # end subsection handling fmcs text drawing

	return(plot_target)
	} # end subfunction: handle_text

handle_segs = function(sdf, plot_target, alpha_edge = 0.5, edge_policy = default_edge_policy(), bond_dist_offset = 0.05, fmcs_result=NULL) {
	# handles drawing of edges on plot_target
	# Args:
	#	sdf: from main function
	#	plot_target: must be a ggplot2 object
	#	alpha_edge: defines transparency of bonds
	#	edge_policy: a policy on which segments to draw, and what color
	#	bond_dist_offset: radial distance to space double or triple bonds by
	#	fmcs_result: if non-null, maximum common substructure in SDF will
	#		be highlighted. If non-null, edge_policy MUST support edge_policy['fmcs']
	
	edge_matrix = cbind(	atomblock(sdf)[bondblock(sdf)[,1], 1:2],
				atomblock(sdf)[bondblock(sdf)[,2], 1:2],
				bondblock(sdf)[,3]
				)

	# evaluate edge pairs to determine if any should be skipped
	node_names = gsub('[^[:alpha:]]', '', rownames(atomblock(sdf)))
	tofrom_matrix = cbind(	node_names[bondblock(sdf)[,1]],
				node_names[bondblock(sdf)[,2]]
				)
	tofrom_matrix = cbind(tofrom_matrix, paste(tofrom_matrix[,1], tofrom_matrix[,2], sep='-'))

	color_vec = cbind(	match(tofrom_matrix[,1], names(edge_policy)),
				match(tofrom_matrix[,2], names(edge_policy)),
				match(tofrom_matrix[,3], names(edge_policy))
				)

	color_vec[is.na(color_vec)] = match('default', names(edge_policy))
	# condense into one vector of max values
	color_vec = apply(color_vec, 1, max)
	color_vec = edge_policy[color_vec]

	if (!is.null(fmcs_result)) {
		fmcs_search_matrix = bondblock(sdf)[,1:2]

		for (i in 1:length(color_vec)) {
			# if both origin atom and target atom are in fmcs
			is_in_fmcs= all(c(	fmcs_search_matrix[i,1],
						fmcs_search_matrix[i,2]) %in% fmcs_result)
			if(is_in_fmcs) {
				# then mark its index in color_vec with fmcs color
				color_vec[i] = edge_policy['fmcs']
				}

			} # end iteration through rows of color_vec using i
		} # end if: fmcs_result is non-null

	# cut out edges that have been marked for skipping
	tofrom_matrix = tofrom_matrix[	color_vec != 'skip',]
	edge_matrix = 	edge_matrix[	color_vec != 'skip',]
	color_vec = 	color_vec[	color_vec != 'skip']

	# evaluate edge weights, split those with str > 1 
	strong_edge_matrix = edge_matrix[edge_matrix[,5] > 1, ]
	strong_color_vec = color_vec[edge_matrix[,5] > 1]
	i = 1; offset_vec = c(1, -1, 1, -1)
	catting_matrix = matrix(data=0, nrow=1, ncol=4)
	catting_color_vec = c(0)
	while (i <= dim(strong_edge_matrix)[1]) {
		dxdy = c(strong_edge_matrix[i,1] - strong_edge_matrix[i,3], strong_edge_matrix[i,2] - strong_edge_matrix[i,4])
		scaled_dxdy = rev((dxdy/sqrt(sum(dxdy^2)))) * bond_dist_offset
		dp1 = strong_edge_matrix[i,1:4] + scaled_dxdy*  1*offset_vec
		dp2 = strong_edge_matrix[i,1:4] + scaled_dxdy* -1*offset_vec
		catting_matrix = rbind(catting_matrix, dp1, dp2)
		catting_color_vec = cbind(catting_color_vec, strong_color_vec[i], strong_color_vec[i])
		
		# if triple+ bond
		if (strong_edge_matrix[i,5] >= 3) {
			# reinsert old bond's xy values
			catting_matrix = rbind(catting_matrix, strong_edge_matrix[i,1:4])
			catting_color_vec = cbind(catting_color_vec, strong_color_vec[i])
			}
		i = i +1
		} # end iteration through rows of strong_edge_matrix using i
	catting_matrix = cbind(catting_matrix, 1)
	if (dim(catting_matrix)[1] > 1) {
		catting_color_vec = catting_color_vec[2:length(catting_color_vec)]
		catting_matrix = catting_matrix[2:dim(catting_matrix)[1],]

		color_vec = color_vec[edge_matrix[,5] <= 1]
		color_vec = c(color_vec, catting_color_vec)
		edge_matrix = edge_matrix[edge_matrix[,5] <= 1,]
		edge_matrix = rbind(edge_matrix, catting_matrix)
		}

	# finalize edge matrix
	rownames(edge_matrix) = NULL
	edge_matrix = data.frame(edge_matrix)
	plot_target = plot_target + geom_segment(alpha = alpha_edge, data=edge_matrix, aes(x=C1, y=C2, xend=C1.1, yend=C2.1), color = color_vec)
	return(plot_target)
	}

draw_sdf = function(	sdf,
			filename = 'test.jpg',
			alpha_edge = 0.5,
			alpha_node = 1.0,
			numbered=FALSE,
			font_size=5,
			node_vertical_offset = 0,
			node_background_color = FALSE,
			bgcolor = rgb(1,1,1,1),
			bgraster = NULL,
			node_policy = default_node_policy(),
			edge_policy = default_edge_policy(),
			bond_dist_offset = 0.05,
			fmcsR_sdf = NULL
			) {
	# draws an SDF
	# Args:
	# 	sdf: must support atomblock() and bondblock(). supports
	#		 length() > 1
	#	alpha_edge: value in range [0=least, 1=most], describes translucence of edge colors
	#	alpha_node: value in range [0=least, 1=most], describes translucence of node colors
	#	numbered: FALSE = node text does not include numbering, TRUE or 1 = numbering style 1,
	#		2 = numbering style 2
	#	font_size: absolute font size. Reduce this value when plotting multiple SDFs
	#	node_vertical_offset: permits text drawing of SDFs to be shifted upward
	#	bgcolor: rgb object that defines base color of background of plot
	#	bgraster: a readPNG object or a path to an object that can be read as PNG. Will
	#		be pasted as first/background layer of image
	#	node_policy: string mapping that defines what nodes to draw, and what color to draw
	#	edge_policy: string mapping that defines what edges to draw, and what color to draw
	#		policies require 'default' to be mapped to a color, ie at minumum set both equal to
	#		c('default'='black')
	#	bond_dist_offset: defines absolute space between double+ bond lines
	#	fmcsR_sdf: an SDF which supports fmcsr(sdf, fmcsR_sdf). Maximum common substructure will
	#		be highlighted in drawing of SDF. If non-null, node_policy and edge_policy MUST
	#		support node_policy['fmcs'] and edge_policy['fmcs']
	#
	#	Sample usage:
	#	> library(ChemmineR)
	#	> data(sdfsample)
	#	> source('draw_sdf.r')
	#	> draw_sdf(sdfsample[[1]], alpha_edge = 0.25, font_size=5)
	#	>
	#	> library(fmcsR)
	#	> draw_sdf(sdfsample[[2]], filename='././file_one.png', bgraster = 'demo_raster.png', fmcsR_sdf = sdfsample[[2]])
	
	if (length(sdf) > 1) {
		return(concatenate_plots(sdf, filename, alpha_edge, alpha_node, numbered, font_size, node_vertical_offset, node_background_color, bgcolor, bgraster, node_policy, edge_policy, fmcsR_sdf))
		# todo add in all args as we develop this function
		} # end if: multiple SDF

	# begin plotting
	
	plot_target = ggplot(data=data.frame(0), aes())
	if (!is.null(bgraster)) {
		plot_target = handle_raster(plot_target, bgraster)
		} # end if: raster was passed

	fmcs_result = NULL
	if (!is.null(fmcsR_sdf)) {
		require(fmcsR)
		fmcs_result = fmcs(sdf, fmcsR_sdf)
		fmcs_result = fmcs_result@mcs1$mcs1$CMP1_fmcs_1
		}
	plot_target = handle_segs(sdf, plot_target, alpha_edge, edge_policy, bond_dist_offset, fmcs_result)
	plot_target = handle_text(sdf, plot_target, alpha_node, numbered, font_size, node_vertical_offset, node_background_color, node_policy, fmcs_result)

	# fix background
	plot_target = plot_target + theme(panel.background = element_rect(fill=bgcolor, color='black'))
	plot_target = plot_target + theme(legend.position='none')
	plot_target = plot_target + theme(axis.title.x = element_blank())
	plot_target = plot_target + theme(axis.title.y = element_blank())

	if (!is.null(filename)) {
		ggsave(filename=filename, plot = plot_target)
		} # end if: filename is nonnull

	return(plot_target)

	} # end function draw_sdf
