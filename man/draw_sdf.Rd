\name{draw_sdf}
\alias{draw_sdf}
\title{
	draw_sdf
}
\description{
	Draws an sdf object in the 2D plane using ggplot2 library. Permits customization of bond colors and atom colors.
}
\usage{
	draw_sdf(sdf, filename = "test.jpg", alpha_edge = 0.5, alpha_node = 1, numbered = FALSE, font_size = 5, node_vertical_offset = 0, node_background_color = FALSE, bgcolor = rgb(1, 1, 1, 1), bgraster = NULL, node_policy = default_node_policy(), edge_policy = default_edge_policy(), bond_dist_offset = 0.05, fmcsR_sdf = NULL)
}
\arguments{
  \item{sdf}{
		An instance of a SDF or list of SDFs
	}
  \item{filename}{
		Filename to save image to. Defaults to 'test.jpg'. If set to NULL, does not save image.
	}
  \item{alpha_edge}{
		alpha of bonds in your image. Defaults to 0.5. 0 is fully transparent, 1 is fully opaque.
	}
  \item{alpha_node}{
		alpha of atoms in your image. Defaults to 1.0.
	}
  \item{numbered}{
		If 1 or TRUE, displays numbering of atoms at their location. If 2, displays a second numbering style.
	}
  \item{font_size}{
		Controls size of text to be displayed at atom locations. Beware when plotting multiple SDFs in one image. Ggplot will still scale fonts as if text is being plotted in one image.
	}
  \item{node_vertical_offset}{
		Upward shift of atom text. Upward shit is in SDF units, not ggplot units.
	}
  \item{bgcolor}{
		An rgb(r,g,b,alpha) or similar object. produces a background of the specified color.
	}
  \item{node_background_color}{
		 A common color as a text string (e.g. 'white', 'pink') or an rgb(r,g,b,alpha).  Draws a filled circle of the color specified before drawing text over each node.
	}
  \item{bgraster}{
		A readPNG object or a path to an object that can be understood using readPNG. Will be used as background.
	}
  \item{node_policy}{
		Mapping that defines how atom strings should be displayed. Simplest would be c('default'='black')
	}
  \item{edge_policy}{
		Mapping that defines how bonds should be displayed. Simplest is c('default'='black'), though this will display all Hydrogen bonds as well.
	}
  \item{bond_dist_offset}{
		Defines space between double or triple bonds, in SDF units.
	}
  \item{fmcsR_sdf}{
		A second SDF object to run fmcsR on. 
	}
}
\details{
	Requires ggplot2. Additional features require grid, gridExtra, fmcsR, or png. Most matrix operations vectorized.
}
\value{
	Returns a ggplot2 object. Calling draw_sdf(...) rather than assigning it will result in R trying to print a ggplot2 object.
}
\author{
	John A. Sharifi
}


\examples{
	library(ChemmineR) 		# if not already imported
	data(sdfsample)
	draw_sdf(sdfsample[[1]])
}
