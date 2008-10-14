# options
if (is.null(getOption('chemminer.max.upload')))
    options(chemminer.max.upload = 100)
    
# supported elements
.elements <- list(
    R=0,
    H=1,
    He=2,
    Li=3,
    Be=4,
    B=5,
    C=6,
    N=7,
    O=8,
    F=9,
    Ne=10,
    Na=11,
    Mg=12,
    Al=13,
    Si=14,
    P=15,
    S=16,
    Cl=17,
    Ar=18,
    K=19,
    Ca=20,
    Sc=21,
    Ti=22,
    V=23,
    Cr=24,
    Mn=25,
    Fe=26,
    Co=27,
    Ni=28,
    Cu=29,
    Zn=30,
    Ga=31,
    Ge=32,
    As=33,
    Se=34,
    Br=35,
    Kr=36,
    Rb=37,
    Sr=38,
    Y=39,
    Zr=40,
    Nb=41,
    Mo=42,
    Tc=43,
    Ru=44,
    Rh=45,
    Pd=46,
    Ag=47,
    Cd=48,
    In=49,
    Sn=50,
    Sb=51,
    Te=52,
    I=53,
    Xe=54,
    Cs=55,
    Ba=56,
    La=57,
    Ce=58,
    Pr=59,
    Nd=60,
    Pm=61,
    Sm=62,
    Eu=63,
    Gd=64,
    Tb=65,
    Dy=66,
    Ho=67,
    Er=68,
    Tm=69,
    Yb=70,
    Lu=71,
    Hf=72,
    Ta=73,
    W=74,
    Re=75,
    Os=76,
    Ir=77,
    Pt=78,
    Au=79,
    Hg=80,
    Tl=81,
    Pb=82,
    Bi=83,
    Po=84,
    At=85,
    Rn=86,
    Fr=87,
    Ra=88,
    Ac=89,
    Th=90,
    Pa=91,
    U=92,
    Np=93,
    Pu=94,
    Am=95,
    Cm=96,
    Bk=97,
    Cf=98,
    Es=99,
    Fm=100,
    Md=101,
    No=102,
    Lr=103,
    Rf=104,
    Db=105,
    Sg=106,
    Bh=107,
    Hs=108,
    Mt=109,
    Ds=110,
    Rg=111)

# read a sdf file and returns the atoms with bond list
# argument: the filename, which can also be a connection
# return: the list containing the vector of atoms and the vector of bonds
# this is an internal function and is not intended to be used by users
.parse <- function(file_name, skip=3)
{
    if ("connection" %in% class(file_name))
        con <- file_name
    else
        con <- file(file_name, open="r")
    num_atoms <- 0
    num_bonds <- 0
    # read header line, skip a few lines (given in `skip')
    if (skip > 1)
        readLines(con=con, n=skip)
    line <- readLines(con=con, n=1)
    num_atoms <- as.integer(substr(line, 1, 3))
    num_bonds <- as.integer(substr(line, 4, 6))

    if (length(num_atoms) == 0 || num_atoms == 0)
        return(list(atoms=vector(), bonds=list(u=vector(), v=vector())))

    # parse atom block
    temp <- readLines(con=con, n=num_atoms)
    atoms <- gsub(' ', '', substr(temp, 32, 34))

    # parse bond block
    # bonds <- array(0, c(num_bonds, 3))
    temp <- readLines(con=con, n=num_bonds)
    u <- as.integer(substr(temp, 1, 3))
    v <- as.integer(substr(temp, 4, 6))
    t <- as.integer(substr(temp, 7, 9))
    bonds <- list(u=u, v=v, t=t)

# return: a compound object:
#   keys:
#       'atoms': a vector of atoms
#       'bonds': a list indexed by 'u', 'v', and 't'. Each entry will
#       be a vector 

    # close connection if it is opend here
    if (! "connection" %in% class(file_name))
        close(con)

    return(list(atoms=atoms, bonds=bonds))
}

# generates distance matrix for every pair of atoms in compound. Distance is
# the length of shortest path connecting them.
# Returns the distance matrix giving distance between any 2 atoms in the
# compound
# this is an internal function and is not intended to be used by users
.all_pairs_dist <- function(cmp)
{
    s <- length(cmp[['atoms']])
    dis <- array(1024, c(s, s))
    # for every edge, init the distance to 1
    for (i in 1:length(cmp[['bonds']][['u']])) {
        dis[cmp[['bonds']][['u']][i], cmp[['bonds']][['v']][i]] <- 1
        dis[cmp[['bonds']][['v']][i], cmp[['bonds']][['u']][i]] <- 1
    }
    for (k in 1:s) {
        for (i in 1:(s-1)) {
        if (i != k)
            for (j in (i + 1):s) {
            if (j != k)
                if (dis[i, k] + dis[k, j] < dis[i, j]) {
                    dis[i, j] <- dis[i, k] + dis[k, j]
                    dis[j, i] <- dis[i, j]
                }
            }
        }
    }

    return(dis)
}

# generates atom descriptor for atom #i. Descriptor is encoded as a number.
# this is an internal function and is not intended to be used by users
# see also: .interpret_atom
.gen_atom_desc <- function(cmp, i)
{
    neighbors <- 0 # neighbours (non hydrogens)
    pi_electrons <- 0 # number of pi electrons
    num_bonds <- length(cmp[['bonds']][['u']])

    # for every bond, update atom properties
    for (k in 1:num_bonds) {
        if (cmp[['bonds']][['u']][k] == i || cmp[['bonds']][['v']][k] == i) {
            if (cmp[['atoms']][ cmp[['bonds']][['u']][k] ] != 'H' &&
                    cmp[['atoms']][ cmp[['bonds']][['v']][k] ] != 'H') {
                neighbors <- neighbors + 1
                # single to triple bonds
                if (cmp[['bonds']][['t']][k] < 4)
                    pi_electrons <- pi_electrons + cmp[['bonds']][['t']][k] - 1
                # arromatic bonds
                else
                    pi_electrons <- pi_electrons + 0.5
            }
        }
    }

    pi_electrons <- floor(pi_electrons)
    atom_desc <- .elements[[cmp[['atoms']][i]]] * 2^6 + neighbors * 2^3 +
            pi_electrons
    return(atom_desc)
}

# generate atom descriptors for each atom of the compound. 
# this is an internal function and is not intended to be used by users
.gen_all_atom_desc <- function(cmp)
{
    num_atoms <- length(cmp[['atoms']])
    desc <- array()
    for (i in 1:num_atoms) {
        desc[[i]] <- .gen_atom_desc(cmp, i)
    }
    return(desc)
}

# generates atom pair descriptor for a compound. It will encode all atom pairs.
# descriptor = atom descriptors + atom pair information
# this is an internal function and is not intended to be used by users
.gen_atom_pair <- function(cmp)
{
    num_bonds <- length(cmp[['bonds']][['u']])
    if (num_bonds == 0)
        return(c())
    dis <- .all_pairs_dist(cmp)
    num_atoms <- length(cmp[['atoms']])
    desc <- c()
    tail <- 1

    a_desc <- .gen_all_atom_desc(cmp)
    # for every pair of atoms
    for (i in 1:num_atoms) {
        for (j in i:num_atoms) {
            if (i != j && dis[i,j] < 128 && cmp[['atoms']][i] != 'H' &&
                    cmp[['atoms']][j] != 'H') {
                # get the atom descriptor for the two atoms
                desc_i <- a_desc[[i]]
                desc_j <- a_desc[[j]]
                # construct descriptor for the atom pair
                if (desc_i >= desc_j)
                    desc[tail] <- desc_i * 2^20 + dis[i,j] * 2^13 + desc_j
                else
                    desc[tail] <- desc_j * 2^20 + dis[i,j] * 2^13 + desc_i
                tail <- tail + 1
            }
        }
    }
    if (length(desc) == 0)
        return(desc)
    else 
        # remove duplicated atom pair descriptors
        return(.factor_to_vector(as.factor(desc)))
}

# similarity model
# input: descriptors for two compounds. both are vectors
# output: similarity
cmp.similarity <- function(a, b, mode=1, worst=0)
{
    # check argument, make sure they are vectors and not lists
    if (class(b) != "numeric" || class(a) != "numeric") {
        stop("Both arguments must be descriptors in the form of vectors.\n",
    "Did you use \"[]\" instead of \"[[]]\" to index the descriptor ",
    "database?")
    }

    if (mode == 1 && worst != 0) {
        # estimate upper bound
        i <- length(a)
        j <- length(b)
        u_bond <- min(i, j) / max(i, j)
        if (u_bond < worst)
            return(0)
    }
    if (mode == 0) {
        # if no mode is given, guess a mode based on compound size
        # difference 
        if (length(a) < .25 * length(b) || length(a) > 4 * length(b))
            mode <- 3
        else
            mode <- 1
    }
    if (mode == 1) 
        return (length(intersect(a, b)) / length(union(a, b)))
    else if (mode == 2)
        return (length(intersect(a, b)) / min(length(a), length(b)))
    else if (mode == 3) {
        s <- length(intersect(a, b)) / min(length(a), length(b))
        return (s^3)
    }
}

# generate (atom pair) descriptor from sdf file encoding ONE compound. If the
# file encodes many compounds, only the first one will be parsed. Use
# cmp.parse/sdfs_to_desc to parse multiple compounds.
# returns descriptor for one compound as vector.
sdf_to_desc <- function(filename)
{
    cmp <- .parse(filename)
    return (.gen_atom_pair(cmp))
}
# alias
cmp.parse1 <- sdf_to_desc

# generate descriptors for all compounds in an sdf file
# returns a list of descriptors of compounds. The list contains two names:
# descdb: the descriptor db, itself being a vector of descriptors.
# names: the names of compounds as a vector.
sdfs_to_desc <- function(filename, quiet=FALSE)
{
    # count compounds first
    con <- file(filename, open="r")
    cat("counting number of compounds in sdf...\n")
    n_compounds <- 0
    while (TRUE) {
        line <- readLines(con, n=1)
        if (length(line) == 0) {
            # reaching the last line
            if (prev_line != "$$$$") {
                n_compounds <- n_compounds + 1
                warning("The SDF file does not end with $$$$!")
            }
            break;
        }
        if (line == "$$$$")
            n_compounds <- n_compounds + 1
        prev_line <- line
    }
    cat(n_compounds, " compounds found\n")
    close(con)

    # real parsing
    con <- file(filename, open="r")
    descdb <- list()
    cids <- array()
    sdf_seg <- array()
    flag <- 1
    id <- 1
    cur_line <- 0
    cmp_start <- cur_line

    # compound id
    cid <- readLines(con, n=1)
    cur_line <- cur_line + 1
    if (length(cid) == 0) cid <- ""

    while (flag == 1) {
        # header
        cmp <- .parse(con, skip=2)
        cur_line <- cur_line + 2 + length(cmp[["atoms"]]) +
            length(cmp[["bonds"]]$u) + 1
        if (length(cmp[["atoms"]]) == 0)
            flag <- 0
        else {
            # skip to the end of cmp entry
            tmp <- readLines(con, n=1)
            cur_line <- cur_line + 1

            while (length(tmp) != 0 && tmp != "$$$$") {
                tmp <- readLines(con, n=1)
                cur_line <- cur_line + 1
            }
        
            # update db
            descdb[[id]] <- .gen_atom_pair(cmp)
            cids[[id]] <- cid
            sdf_seg[[id]] <- cmp_start
            if (! quiet) {
                cid_p <- substr(cid, 1, 20)
                if (cid_p == cid) 
                    msg <- paste("\r", id, "/", n_compounds, " parsed (", 
                        cid_p, ")",
                        " now at line ", cur_line,
                        '                                                     ',
                        sep='')
                else
                    msg <- paste("\r", id, "/", n_compounds, " parsed (", 
                        cid_p, "...)",
                        " just passed line ", cur_line,
                        '                                                     ',
                        sep='')
                cat(substr(msg, 1, 79))
            }
            id <- id + 1

            # proceed to next compound, which might not exist
            cmp_start <- cur_line
            cid <- readLines(con, n=1)
            cur_line <- cur_line + 1
            if (length(cid) == 0)
                flag <- 0
        }

    }
    cat(
    "\nyou can use save(..., file='...', compress=TRUE) to save the database\n"
    )
    close(con)
    return(list(descdb=descdb, cids=cids, sdfsegs=sdf_seg, source=filename))
}
# alias
cmp.parse <- sdfs_to_desc

# search db for similar compounds. `db' gives the database returned by
# `cmp.parse' or `sdfs_to_desc'. `query' gives the descriptor of one compound,
# either returned by `cmp.parse1/sdf_to_desc' or a reference to one instance in
# db such as db$descdb[[123]] two types of cutoff: score cutoff (<=1) or count
# cutoff (>1)
cmp.search <- function(db, query, cutoff=0.5, return.score=FALSE, quiet=FALSE,
        mode=1, visualize=FALSE, visualize.browse=TRUE, visualize.query=NULL)
{
    scores <- array()
    ids <- array()
    pos <- 1    # tail of entries
    size <- length(db$cids)
    perstep <- round(size / 25)
    for (i in 1:size) {
        if (!quiet && i %% perstep == 0) {
            steps <- i / perstep
            .progress_bar(paste(min(steps*4, 100), "%", collapse=""))
        }
        if (cutoff <= 1)
            score <- cmp.similarity(db$descdb[[i]], query, mode=mode,
                worst=cutoff)
        else if (pos - 1 < cutoff)
            score <- cmp.similarity(db$descdb[[i]], query, mode=mode, worst=0)
        else
            score <- cmp.similarity(db$descdb[[i]], query, mode=mode,
                    worst=scores[[pos - 1]])
        if (cutoff <= 1) {
            if (score >= cutoff) {
                scores[[pos]] <- score
                ids[[pos]] <- i
                pos <- pos + 1
            }
        } else {
            # using naive insertion sort
            if (pos - 1 < cutoff)
                pos <- pos + 1
            else if (scores[[pos - 1]] >= score)
                next

            pos_ <- pos - 2
            while (pos_ > 0 && scores[[pos_]] < score) {
                scores[[pos_ + 1]] <- scores[[pos_]]
                ids[[pos_ + 1]] <- ids[[pos_]]
                pos_ <- pos_ - 1
            }
            pos_ <- pos_ + 1
            scores[[pos_]] <- score
            ids[[pos_]] <- i
        }
    }
    if (!quiet) cat("\n")
    if (cutoff <= 1) {
        order <- sort(scores, index=TRUE)[["ix"]]
        order <- order[length(order):1]
    } else
        order <- 1:length(scores)

    if (visualize) {
        notes <- as.character(scores[order])
        names(notes) <- rep("similarity", length(notes))
        if (is.null(visualize.query))
            url <- sdf.visualize(db, ids[order], extra=notes,
                browse=visualize.browse)
        else
            url <- sdf.visualize(db, ids[order], reference.sdf=visualize.query,
                extra=notes, browse=visualize.browse)
        if (!visualize.browse)
            print(url)
    }
    if (return.score)
        return(data.frame(ids=ids[order], scores=scores[order]))
    else
        return(ids[order])
}

# view sdfs in ChemMine
.sdf.visualize.max.files = 100
sdf.visualize <- function(db, cmps, extra=NULL, reference.sdf=NULL,
    reference.note=NULL, browse=TRUE, quiet=TRUE)
{
    if (is.null(getOption('chemminer.max.upload')))
        max.upload = 100
    else
        max.upload = getOption('chemminer.max.upload')
    if (length(cmps) > max.upload)
        stop(paste('\n\tYou can visualize at most ', 
            max.upload,
            ' compounds.\n',
            '\tYou supplied ', length(cmps), ' compounds.\n',
            '\tSending too many compounds will take a long time and much\n',
            '\tresource on the server, and may also crash your browser.\n',
            '\tYou can reset the limit by setting the `chemminer.max.upload` option.\n',
            sep=''))

    # read the reference file if there is any
    if (!is.null(reference.sdf))
        reference.sdf <- .read.one.sdf(reference.sdf)

    if (!is.null(extra) && length(extra) != length(cmps))
        stop(
        '\n\tthe indices and the extra information have different lengths.\n',
        '\tYou supplied ', length(cmps), ' compounds.\n',
        '\t`extra\' has a length of ', length(extra)
        )

    sdfs <- sdf.subset(db, cmps)
    cids <- db$cids[cmps]

    # build query
    query <- list(sdf=sdfs, cids=paste(cids, collapse="\1"))
    if (! is.null(extra)) {
        # if an entry in extra is a data frame, format it
        .extra <- c()
        for (i in extra) {
            if (class(i) == 'data.frame' || class(i) == 'matrix')
                .extra <- c(.extra, .data.frame.to.str(i))
            else
                .extra <- c(.extra, paste(as.character(i)))
        }
        query <- c(query, list(extra=paste(.extra, collapse="\1")))
        if (! is.null(names(extra)))
            query <- c(query, list(names=paste(names(extra), collapse="\1")))
    }
    if (! is.null(reference.sdf)) {
        query <- c(query, list(referencesdf=reference.sdf))
        if (! is.null(reference.note)) {
            if (class(reference.note) == 'list' && length(reference.note) == 1) 
            {
                    if (! is.null(names(reference.note)))
                        query <- c(query,
                            list(referencenotename=names(reference.note)))
                    reference.note <- reference.note[[1]]
            }
            if (class(reference.note) == 'data.frame' ||
                class(reference.note) == 'matrix')
                query <- c(query, 
                    list(referencenote=.data.frame.to.str(reference.note)))
            else
                query <- c(query, 
                    list(referencenote=paste(as.character(reference.note))))
        }
    }

    # sending
    if (! quiet) cat('sending SDF to ChemMine\n')
    response <- .postToHost("bioweb.ucr.edu",
        "/ChemMineV2/chemminer/postsdfs", query)
    # reading response
    ref <- gsub(
    '.*\\.\\.\\.\\.\\.\\.\\.\\.\\.\\.(.*)\\.\\.\\.\\.\\.\\.\\.\\.\\.\\..*',
    '\\1', response)
    if (! quiet)
        cat('starting your browser (use `options(browser="...")` to', 
        'customize browser)\n')
    url <- paste(
        c("http://bioweb.ucr.edu/ChemMineV2/chemminer/viewsdfs?", ref),
        collapse='')
    if (browse) browseURL(url)
    return(url)
}

# segment SDFs. given indecies of selected compounds, return the concatenation
# of SDFs for them
sdf.subset <- function(db, cmps)
{
    if (!grep('://', db$source) && !file.exists(db$source)) {
        stop("cannot find the source SDF file for this database\n")
    }
    output <- ""
    for (i in cmps)
    {
        sdf <- .extract_single_sdf(db, i)
        output <- paste(output, sdf, sep="")
    }

    return(output)
}

db.subset <- function(db, cmps)
{
    return(list(descdb=db$descdb[cmps], cids=db$cids[cmps],
                sdfsegs=db$sdfsegs[cmps], source=db$source))
}

# helper function for subset_sdf. This function takes an index, and return the
# SDF correponding compound at that index
.extract_single_sdf <- function(db, index)
{
    if (index > length(db$sdfsegs))
        stop("Index out of range in .extract_single_sdf")

    skip <- db$sdfsegs[[index]]
    con <- file(db$source, open="r")
    if (skip > 0) 
        #scan(con, "raw", skip=skip-1, nlines=1, quiet=TRUE)
        readLines(con, skip)
    sdf <- .read.one.sdf(con, collapse=FALSE)

    # replace line 1 with compound index if no name is available
    if (length(grep("[a-zA-z0-9]", sdf[[1]])) == 0)
        sdf[[1]] <- paste("ChemmineR_Unnamed_Compound_", index, sep="")

    # close connection
    close(con)
    return(paste(sdf, sep="", collapse="\n"))
}

# explain an atom pair descriptor
db.explain <- function(desc)
{
    if (length(desc) == 1)
        return(.explain(desc))
    else {
        result <- c()
        for (i in 1:length(desc))
            result <- c(result, .explain(desc[[i]]))
        return(result)
    }
}

.explain <- function(desc)
{
    # numbers of appearance
    occurence <- desc %% 2^7 + 1
    desc <- desc %/% 2^7
    # extract left atom, distance, and right atom
    left <- desc %/% 2^20
    desc <- desc %% 2^20
    dist <- desc %/% 2^13
    desc <- desc %% 2^13
    right <- desc
    # atom properties
    left <- .interpret_atom(left)
    right <- .interpret_atom(right)

    # output
    return (paste(left$symbol, " [", left$neighbors, " neighbor(s),",
                left$pi_electrons, " pi electrons]<---",
                dist,
                "--->", right$symbol, " [", right$neighbors,
                " neighbor(s),", right$pi_electrons, " pi electrons]",
                collapse="", sep=""))
}

# interpret an atom descriptor
# see also: .gen_atom_desc
.interpret_atom <- function(atom_desc)
{
    # atom
    atom_symbol <- names(.elements[.elements == atom_desc %/% 2^6])
    # neighbors
    neighbors <- atom_desc %% 2^6 %/% 2^3
    # pi electron
    pi_electrons <- atom_desc %% 2^3
    return(list(symbol=atom_symbol, neighbors=neighbors,
                pi_electrons=pi_electrons))
}

# convert a factor to a vector to simplify Union and Intersection operations.
# It basically eliminate repeated names by adding suffix. For example, if `a'
# appears twice in the factor, it will be converted to `a1' and `a2'. 
# this is an internal function and is not intended to be used by users
.factor_to_vector <- function(fac)
{
    vec <- vector(mode='numeric')
    index <- 1
    freq <- table(fac)
    keys <- names(freq)
    for (i in 1:length(keys)) {
        key <- keys[[i]]
        if (freq[[key]] == 1) {
            vec[[index]] <- as.numeric(key) * 2^7
            index <- index + 1
        } else {
            for (j in 1:freq[[key]]) {
                vec[[index]] <- as.numeric(key) * 2^7 + j - 1
                index <- index + 1
            }
        }
    }

    return(vec)
}

.progress_bar_int_cnt <- 0
.progress_bar <- function(label=NULL) {
    progress <- c("|", "/", "-", "\\")
    unlockBinding(".progress_bar_int_cnt", environment(.progress_bar))
    .progress_bar_int_cnt <<- .progress_bar_int_cnt %% 4 + 1
    if (is.null(label)) {
        cat("\r", progress[[.progress_bar_int_cnt]])
    } else {
        cat("\r                                           \r")
        cat(progress[[.progress_bar_int_cnt]], label)
    }
}

.postToHost <- function(host, path, data.to.send, port=80)
{
    if(missing(path)) path <- "/"
    ua <- "User-Agent: Mozilla/5.0 (X11; U; Linux i686; en-US; rv:1.8.1.1)"


    boundary <- paste(
            as.integer(runif(20, min=1, max=10)), collapse=""
    )

    header <- NULL
    header <- c(header, paste("POST ", path, " HTTP/1.1", sep=""))
    header <- c(header, paste("Host: ", host, sep=""))
    header <- c(header, "Connection: close")
    header <- c(header, ua)
    header <- c(header, "Accept: */*")
    header <- c(header,
            paste("Content-type: multipart/form-data; boundary=",
                boundary, sep="")
    )
    boundary <- paste('--', boundary, sep='')

    mcontent <- NULL # keeps the content.
    size <- 0

    for (i in 1:length(data.to.send)) {
        value <- data.to.send[[i]]
        key <- names(data.to.send)[i]

        content <- paste(
            boundary, "\n",
            "Content-Disposition: form-data; name=\"", as.character(key), "\"",
            "\n",
            "\n",
            as.character(value), "\n",
            sep=""
        )
        mcontent <- c(mcontent, charToRaw(content))
    }

    size <- length(mcontent) + length(strsplit(boundary, "")[[1]]) + 4
    header <- c(header, paste("Content-length: ", size, "\n\n", sep=""))

    postdata <- c(
        charToRaw(paste(header, collapse="\n")),
        mcontent,
        charToRaw(paste(boundary,"--\n\n",sep=""))
    )

    rm(header,mcontent)

    scon <- socketConnection(host=host,port=port,open="a+b",blocking=TRUE)
    writeBin(postdata, scon, size=1)

    output <- character(0)
    repeat{
        ss <- rawToChar(readBin(scon, "raw", 2048))
        output <- paste(output,ss,sep="")
        if(regexpr("\r\n0\r\n\r\n",ss)>-1) break()
        if(ss == "") break()
    }
    close(scon)
    return(output)
}

.read.one.sdf <- function(filename, collapse=TRUE)
{
    if (length(filename) != 1)
        stop('reference sdf must be a vector of only one entry!')

    if (!("connection" %in% class(filename)) &&
        !("character" %in% class(filename)))
        stop('reference sdf must be an SDF file or one character string of sdf')

    if ("connection" %in% class(filename))
        con <- filename
    # is filename really a filename a the SDF content? check existence of \n
    else if (length(grep("\n", filename)) == 0 &&
        length(grep("\r", filename)) == 0)
        con <- file(filename, open="r")
    else
        con <- textConnection(filename)

    result <- c()
    while (TRUE) {
        line <- readLines(con=con, n=1)
        if (length(line) == 0 || "$$$$" %in% line) {
            if (length(result) == 0)
                stop("empty SDF!")
            result <- c(result, '$$$$\n')
        if (collapse)
            return (paste(result, collapse="\n"))
        else
            return (result)
        } else
            result <- c(result, line)
    }
}

.data.frame.to.str <- function(x)
{
    string <- ''    # have to define this; otherwise will complain it undefined
    con <- textConnection("string", "w")
    write.table(format(x), con, row.names=F, quote=F)
    close(con)
    return(paste(string, collapse='\n'))
} 
