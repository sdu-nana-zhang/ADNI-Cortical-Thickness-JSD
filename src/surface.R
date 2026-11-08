.surface_kernel_cache <- new.env(parent=emptyenv())
.surface_normal_mass <-
function (a, b) 
{
    stopifnot(identical(dim(a), dim(b)), length(a) == length(b), all(a <= b))
    ans <- pnorm(b) - pnorm(a)
    positive <- a > 0
    ans[positive] <- pnorm(a[positive], lower.tail = FALSE) - pnorm(b[positive], lower.tail = FALSE)
    ans
}

.surface_read_scan_original <-
function (path) 
{
    expected <- surface_roi_names()
    rois <- setNames(vector("list", 62L), expected)
    qc <- list()
    sources <- character()
    geom <- list()
    prov <- surface_provenance(path)
    for (hemi in c("lh", "rh")) {
        p <- function(folder, suffix) file.path(path, folder, paste0(hemi, ".", suffix))
        mesh <- surface_read_mesh(p("surf", "white"))
        nv <- nrow(mesh$vertices)
        annot <- surface_read_annot(p("label", "aparc.DKTatlas.annot"))
        cortex <- surface_read_cortex(p("label", "cortex.label"), nv)
        if (length(annot$labels) != nv) 
            stop("Mesh/annotation vertex mismatch")
        names31 <- sub("^[lr]h_", "", expected[startsWith(expected, paste0(hemi, "_"))])
        if (!all(names31 %in% annot$color_table$Name)) 
            stop("DKT31 target labels missing from color table")
        unexpected <- setdiff(unique(annot$labels), c(names31, "unknown", "corpuscallosum", "medialwall"))
        if (length(unexpected)) 
            stop(paste("Unexpected non-DKT31 labels:", paste(unexpected, collapse = ",")))
        area <- surface_vertex_areas(mesh$vertices, mesh$faces)
        x <- setNames(lapply(c("thickness", "sulc", "curv"), function(m) {
            fn <- p("surf", m)
            if (!file.exists(fn)) 
                return(rep(NA_real_, nv))
            a <- surface_read_morph(fn)
            if (length(a) != nv || attr(a, "faces") != nrow(mesh$faces)) 
                stop("Morphometry/mesh dimensions mismatch")
            as.numeric(a)
        }), c("thickness", "sulc", "curv"))
        if (!prov$WhiteCurvatureCommandVerified) 
            x$curv[] <- NA_real_
        available <- vapply(x, function(v) any(is.finite(v)), logical(1))
        valid <- cortex & is.finite(area) & area > 0
        for (m in names(x)[available]) valid <- valid & is.finite(x[[m]]) & if (m == "thickness") 
            x[[m]] > 0
        else TRUE
        for (nm in names31) {
            r <- paste0(hemi, "_", nm)
            allidx <- annot$labels == nm & cortex
            ii <- which(allidx & valid)
            totalarea <- sum(area[allidx & is.finite(area) & area > 0])
            validarea <- sum(area[ii])
            frac <- if (totalarea > 0) 
                validarea/totalarea
            else 0
            neff <- if (length(ii)) 
                sum(area[ii])^2/sum(area[ii]^2)
            else 0
            pass <- length(ii) >= 30L && neff >= 20 && frac >= 0.9
            reason <- if (pass) 
                ""
            else paste(c(if (length(ii) < 30) "fewer_than_30_vertices", if (neff < 20) "area_effective_n_below_20", 
                if (frac < 0.9) "retained_area_below_90pct"), collapse = ";")
            rois[[r]] <- list(x = lapply(x, function(v) v[ii]), area = as.numeric(area[ii]), vertex_index_zero = ii - 
                1L, metric_available = available, technical_pass = pass, reason = reason)
            qc[[r]] <- data.frame(ROI = r, N_cortex = sum(allidx), N_retained = length(ii), Area_total = totalarea, 
                Area_retained = validarea, RetainedAreaFraction = frac, AreaEffectiveN = neff, InvalidThickness = sum(allidx & 
                  (!is.finite(x$thickness) | x$thickness <= 0)), InvalidSulc = sum(allidx & !is.finite(x$sulc)), 
                InvalidCurv = sum(allidx & !is.finite(x$curv)), ThicknessAtLeast5 = sum(allidx & is.finite(x$thickness) & 
                  x$thickness >= 5), TechnicalPass = pass, Reason = reason)
        }
        sources <- c(sources, p("surf", "white"), p("surf", "pial"), p("surf", "thickness"), p("surf", 
            "sulc"), p("surf", "curv"), p("label", "aparc.DKTatlas.annot"), p("label", "cortex.label"))
        geom[[hemi]] <- list(N_vertices = nv, N_faces = nrow(mesh$faces), SurfaceArea = sum(area), ZeroAreaFaces = attr(area, 
            "zero_area_faces"))
    }
    sources <- c(sources, file.path(path, "scripts", c("recon-all.log", "recon-surf.log", "recon-surf.done", 
        "recon-all.done", "build-stamp.txt")))
    sources <- sources[file.exists(sources)]
    info <- file.info(sources)
    manifest <- data.frame(RelativePath = substring(sources, nchar(path) + 2L), Bytes = info$size, Modified = as.character(info$mtime), 
        MD5 = unname(tools::md5sum(sources)))
    list(ScanKey = basename(path), SourcePath = path, ROI = rois, QC = do.call(rbind, qc), Geometry = geom, 
        Provenance = prov, Manifest = manifest, ExtractionVersion = "surface_dkt62_area_commonvalid_v1", 
        ExtractedAt = format(Sys.time(), tz = "UTC", usetz = TRUE))
}

adni_empty_volumes <-
function () 
{
    d <- adni_volume_dictionary()
    nm <- c(paste0(c(d$globals, d$structures, d$pairs, "CerebralGrayDerived", "SixVentriclesSum"), "_stats_mm3"), 
        paste0("DKT62Volume_stats_mm3__", surface_roi_names()), paste0("DKT62Area_stats_mm2__", surface_roi_names()), 
        paste0("Volume_", adni_label_ids()))
    setNames(rep(NA_real_, length(nm)), nm)
}

adni_error <-
function (e) 
conditionMessage(e)

adni_extract_volumes <-
function (path, read_label_volumes = TRUE) 
{
    d <- adni_volume_dictionary()
    values <- adni_empty_volumes()
    errors <- character()
    long <- list()
    state <- list(AsegStatsParsed = FALSE, DKTLeftStatsParsed = FALSE, DKTRightStatsParsed = FALSE, AsegStatsSource = "", 
        DKTLeftStatsSource = "", DKTRightStatsSource = "", GlobalVolumesComplete = FALSE, Core14VolumesComplete = FALSE, 
        DKT62VolumesComplete = FALSE, LabelVolumesAvailable = FALSE, LabelVolumeSource = "", LabelVoxelVolume_mm3 = NA_real_, 
        LabelSourceMatchesOriginal1996 = FALSE, StatsVolumeError = "", StatsFastSurferVersion = NA_character_)
    aseg <- file.path(path, "stats/aseg.stats")
    if (file.exists(aseg)) 
        tryCatch({
            a <- adni_stats_table(aseg)
            tab <- a$table
            if (!all(c("SegId", "Volume_mm3") %in% names(tab)) || anyDuplicated(tab$SegId)) 
                stop("aseg分区字段缺失或SegId重复。")
            measures <- strsplit(sub("^# Measure ", "", a$lines[grepl("^# Measure ", a$lines)]), ",", 
                fixed = TRUE)
            for (nm in d$globals) {
                q <- measures[vapply(measures, function(x) length(x) >= 5L && trimws(x[2]) == nm, logical(1))]
                if (length(q) > 1L) 
                  stop("全局体积重复：", nm)
                if (length(q)) {
                  if (!trimws(q[[1]][5]) %in% c("mm^3", "mm3")) 
                    stop("全局体积单位异常：", nm)
                  v <- as.numeric(trimws(q[[1]][4]))
                  if (!is.finite(v) || v < 0) 
                    stop("全局体积无效：", nm)
                  values[paste0(nm, "_stats_mm3")] <- v
                }
            }
            for (i in seq_along(d$ids)) {
                j <- match(d$ids[i], tab$SegId)
                if (!is.na(j)) {
                  v <- as.numeric(tab$Volume_mm3[j])
                  if (!is.finite(v) || v < 0) 
                    stop("结构体积无效。")
                  values[paste0(d$structures[i], "_stats_mm3")] <- v
                }
            }
            add <- function(names) {
                a <- values[paste0(names, "_stats_mm3")]
                if (any(!is.finite(a))) 
                  NA_real_
                else sum(a)
            }
            for (i in seq_along(d$pairs)) values[paste0(d$pairs[i], "_stats_mm3")] <- add(c(d$left[i], 
                d$right[i]))
            values["CerebralGrayDerived_stats_mm3"] <- add(c("CortexVol", "SubCortGrayVol"))
            values["SixVentriclesSum_stats_mm3"] <- add(c("LeftLateralVentricle", "RightLateralVentricle", 
                "LeftInferiorLateralVentricle", "RightInferiorLateralVentricle", "ThirdVentricle", "FourthVentricle"))
            state$AsegStatsParsed <- TRUE
            state$AsegStatsSource <- "stats/aseg.stats"
            ver <- grep("^# FastSurfer_version ", a$lines, value = TRUE)
            if (length(ver)) 
                state$StatsFastSurferVersion <- trimws(sub("^# FastSurfer_version ", "", tail(ver, 1L)))
            long[[1]] <- data.frame(Source = "aseg.stats", tab, check.names = FALSE)
        }, error = function(e) errors <<- c(errors, paste("aseg:", adni_error(e))))
    else errors <- c(errors, "缺少stats/aseg.stats")
    for (hemi in c("lh", "rh")) {
        candidates <- paste0("stats/", hemi, c(".aparc.DKTatlas.mapped.stats", ".aparc.DKTatlas.stats"))
        present <- candidates[file.exists(file.path(path, candidates))]
        cap <- if (hemi == "lh") 
            "Left"
        else "Right"
        if (!length(present)) {
            errors <- c(errors, paste("缺少", hemi, "DKT皮层stats"))
            next
        }
        chosen <- present[1L]
        tryCatch({
            z <- adni_stats_table(file.path(path, chosen))
            tab <- z$table
            roi <- sub("^[lr]h_", "", surface_roi_names()[startsWith(surface_roi_names(), paste0(hemi, 
                "_"))])
            if (!all(c("GrayVol", "SurfArea") %in% names(tab))) 
                stop("DKT stats缺少GrayVol或SurfArea。")
            idx <- match(roi, tab$StructName)
            for (i in seq_along(roi)) if (!is.na(idx[i])) {
                v <- as.numeric(tab$GrayVol[idx[i]])
                a <- as.numeric(tab$SurfArea[idx[i]])
                if (any(!is.finite(c(v, a))) || any(c(v, a) < 0)) 
                  stop("DKT皮层体积/面积无效。")
                values[paste0("DKT62Volume_stats_mm3__", hemi, "_", roi[i])] <- v
                values[paste0("DKT62Area_stats_mm2__", hemi, "_", roi[i])] <- a
            }
            state[[paste0("DKT", cap, "StatsParsed")]] <- TRUE
            state[[paste0("DKT", cap, "StatsSource")]] <- chosen
        }, error = function(e) errors <<- c(errors, paste(hemi, adni_error(e))))
    }
    if (isTRUE(read_label_volumes)) {
        candidates <- c("mri/aparc.DKTatlas+aseg.mapped.mgz", "mri/aparc.DKTatlas+aseg.mgz")
        present <- candidates[file.exists(file.path(path, candidates))]
        if (length(present)) 
            tryCatch({
                v <- adni_read_label_mgz(file.path(path, present[1L]))
                values[names(v$volume)] <- v$volume
                state$LabelVolumesAvailable <- TRUE
                state$LabelVolumeSource <- present[1L]
                state$LabelVoxelVolume_mm3 <- v$voxel_volume
                state$LabelSourceMatchesOriginal1996 <- identical(present[1L], candidates[1L])
            }, error = function(e) errors <<- c(errors, paste("标签体积:", adni_error(e))))
        else errors <- c(errors, "缺少DKT+aseg标签MGZ；保留stats体积")
    }
    required <- paste0(c("CortexVol", "SubCortGrayVol", "TotalGrayVol", "eTIV", "LeftHippocampus", "RightHippocampus"), 
        "_stats_mm3")
    state$GlobalVolumesComplete <- state$AsegStatsParsed && all(is.finite(values[required]) & values[required] > 
        0)
    x <- values[paste0(d$structures[1:14], "_stats_mm3")]
    state$Core14VolumesComplete <- state$AsegStatsParsed && all(is.finite(x) & x > 0)
    x <- values[paste0("DKT62Volume_stats_mm3__", surface_roi_names())]
    state$DKT62VolumesComplete <- all(is.finite(x) & x > 0) && state$DKTLeftStatsParsed && state$DKTRightStatsParsed
    state$StatsVolumeError <- paste(errors, collapse = " | ")
    list(values = values, state = state)
}

adni_label_ids <-
function () 
{
    cortex <- c(1002, 1003, 1005:1031, 1034, 1035)
    stopifnot(length(cortex) == 31L)
    c(adni_volume_dictionary()$ids[1:14], cortex, cortex + 1000L)
}

adni_read_label_mgz <-
function (path) 
{
    con <- if (grepl("[.]mgz$", path, ignore.case = TRUE)) 
        gzfile(path, "rb")
    else file(path, "rb")
    on.exit(close(con))
    withCallingHandlers({
        h <- readBin(con, integer(), 7L, size = 4L, endian = "big")
        if (length(h) != 7L || anyNA(h) || h[1] != 1L || any(h[2:4] <= 0L) || h[5] != 1L) 
            stop("MGZ头部异常或不是单帧。")
        n <- prod(as.double(h[2:4]))
        if (n > 1e+08) 
            stop("MGZ体素数过大，需核对头部。")
        ras <- readBin(con, integer(), 1L, size = 2L, endian = "big")
        if (length(ras) != 1L || ras != 1L) 
            stop("MGZ缺少明确空间信息，未假定体素体积为1。")
        delta <- readBin(con, numeric(), 3L, size = 4L, endian = "big")
        directions <- readBin(con, numeric(), 9L, size = 4L, endian = "big")
        center <- readBin(con, numeric(), 3L, size = 4L, endian = "big")
        padding <- readBin(con, raw(), 194L)
        if (length(padding) != 194L || length(delta) != 3L || length(directions) != 9L || length(center) != 
            3L || any(!is.finite(c(delta, directions, center))) || any(delta <= 0)) 
            stop("MGZ空间头部不完整。")
        vv <- abs(det(matrix(directions, 3L) %*% diag(delta)))
        if (!is.finite(vv) || vv <= 1e-10) 
            stop("MGZ体素体积无效。")
        type <- h[6]
        if (!type %in% c(0L, 1L, 3L, 4L)) 
            stop("不支持MGZ数据类型。")
        size <- switch(as.character(type), `0` = 1L, `1` = 4L, `3` = 4L, `4` = 2L)
        v <- readBin(con, if (type == 3L) 
            numeric()
        else integer(), n = n, size = size, signed = type != 0L, endian = "big")
        if (length(v) != n || any(!is.finite(v)) || any(v < 0) || any(v != floor(v)) || max(v) > 1e+05) 
            stop("分区标签被截断或不是合理整数标签。")
        repeat {
            if (!length(readBin(con, raw(), 65536L))) 
                break
        }
        ids <- adni_label_ids()
        counts <- tabulate(as.integer(v) + 1L, nbins = max(c(ids, v)) + 1L)
        list(volume = setNames(counts[ids + 1L] * vv, paste0("Volume_", ids)), voxel_volume = vv, shape = h[2:4])
    }, warning = function(w) stop("MGZ读取警告：", conditionMessage(w), call. = FALSE))
}

adni_stats_table <-
function (path) 
{
    lines <- readLines(path, warn = FALSE)
    head <- grep("^# ColHeaders ", lines, value = TRUE)
    if (length(head) != 1L) 
        stop("stats文件缺少唯一ColHeaders：", path)
    cols <- strsplit(trimws(sub("^# ColHeaders ", "", head)), "[[:space:]]+")[[1]]
    records <- lines[nzchar(trimws(lines)) & !grepl("^[[:space:]]*#", lines)]
    if (!length(records)) 
        stop("stats分区表没有记录：", path)
    tab <- read.table(text = paste(records, collapse = "\n"), header = FALSE, col.names = cols, stringsAsFactors = FALSE, 
        check.names = FALSE, comment.char = "", quote = "")
    if (!"StructName" %in% names(tab) || anyDuplicated(tab$StructName)) 
        stop("stats缺少StructName或结构名称重复：", path)
    list(table = tab, lines = lines)
}

adni_volume_dictionary <-
function () 
{
    list(globals = c("BrainSegVol", "BrainSegVolNotVent", "lhCortexVol", "rhCortexVol", "CortexVol", 
        "SubCortGrayVol", "TotalGrayVol", "eTIV", "VentricleChoroidVol"), ids = c(10, 11, 12, 13, 17, 
        18, 26, 49, 50, 51, 52, 53, 54, 58, 4, 5, 43, 44, 14, 15, 72, 8, 47), structures = c("LeftThalamus", 
        "LeftCaudate", "LeftPutamen", "LeftPallidum", "LeftHippocampus", "LeftAmygdala", "LeftAccumbens", 
        "RightThalamus", "RightCaudate", "RightPutamen", "RightPallidum", "RightHippocampus", "RightAmygdala", 
        "RightAccumbens", "LeftLateralVentricle", "LeftInferiorLateralVentricle", "RightLateralVentricle", 
        "RightInferiorLateralVentricle", "ThirdVentricle", "FourthVentricle", "FifthVentricle", "LeftCerebellumCortex", 
        "RightCerebellumCortex"), pairs = c("HippocampusTotal", "AmygdalaTotal", "ThalamusTotal", "CaudateTotal", 
        "PutamenTotal", "PallidumTotal", "AccumbensTotal", "LateralVentriclesTotal"), left = c("LeftHippocampus", 
        "LeftAmygdala", "LeftThalamus", "LeftCaudate", "LeftPutamen", "LeftPallidum", "LeftAccumbens", 
        "LeftLateralVentricle"), right = c("RightHippocampus", "RightAmygdala", "RightThalamus", "RightCaudate", 
        "RightPutamen", "RightPallidum", "RightAccumbens", "RightLateralVentricle"))
}

surface_extract_scan <-
function (path) 
{
    rel <- c(unlist(lapply(c("lh", "rh"), function(h) c(paste0("surf/", h, ".", c("white", "pial", "thickness", 
        "sulc", "curv")), paste0("label/", h, ".", c("aparc.DKTatlas.annot", "cortex.label"))))), paste0("scripts/", 
        c("recon-all.log", "recon-surf.log", "recon-surf.done", "recon-all.done", "build-stamp.txt")))
    rel <- rel[file.exists(file.path(path, rel))]
    before <- unname(tools::md5sum(file.path(path, rel)))
    names(before) <- rel
    result <- .surface_read_scan_original(path)
    stopifnot(all(result$Manifest$RelativePath %in% rel), all(result$Manifest$MD5 == before[result$Manifest$RelativePath]))
    result$ReadBeforeAfterHashesAgree <- TRUE
    result
}

surface_fast_masses <-
function (rois, metric, config, spacing_divisor = 32L) 
{
    cc <- config[[metric]]
    h <- cc$bandwidth
    delta <- h/spacing_divisor
    indices <- weights <- groups <- list()
    valid <- rep(FALSE, length(rois))
    for (i in seq_along(rois)) {
        a <- rois[[i]]
        if (is.null(a) || !a$technical_pass || !isTRUE(a$metric_available[[metric]])) 
            next
        x <- a$x[[metric]]
        w <- a$area/sum(a$area)
        if (any(!is.finite(x)) || any(!is.finite(w))) 
            stop("Unexpected invalid raw cache")
        u <- x/delta
        lo <- floor(u)
        fraction <- u - lo
        indices[[i]] <- c(lo, lo + 1L)
        weights[[i]] <- c(w * (1 - fraction), w * fraction)
        groups[[i]] <- rep(i, 2L * length(x))
        valid[i] <- TRUE
    }
    if (!any(valid)) 
        return(matrix(NA_real_, length(rois), length(cc$breaks) - 1L, dimnames = list(names(rois), NULL)))
    ix <- unlist(indices, use.names = FALSE)
    ww <- unlist(weights, use.names = FALSE)
    gg <- unlist(groups, use.names = FALSE)
    klo <- min(c(ix, cc$center_index_limits[1] * spacing_divisor/32))
    khi <- max(c(ix, cc$center_index_limits[2] * spacing_divisor/32))
    klo <- floor(klo)
    khi <- ceiling(khi)
    tag <- paste(metric, format(h, digits = 16), spacing_divisor, length(cc$breaks), format(cc$breaks[2], 
        digits = 16), format(cc$breaks[length(cc$breaks) - 1], digits = 16), sep = "|")
    old <- .surface_kernel_cache[[tag]]
    if (!is.null(old) && !identical(old$breaks, cc$breaks)) 
        old <- NULL
    if (!is.null(old) && old$klo <= klo && old$khi >= khi) {
        klo <- old$klo
        khi <- old$khi
        kernel <- old$kernel
    }
    else {
        if (!is.null(old)) {
            klo <- min(klo, old$klo)
            khi <- max(khi, old$khi)
        }
        if (khi - klo > 50000) 
            stop("Extreme data range requires technical review, no clipping performed")
        centers <- (klo:khi) * delta
        lower <- head(cc$breaks, -1)
        upper <- tail(cc$breaks, -1)
        a <- outer(centers, lower, function(x, b) (b - x)/h)
        b <- outer(centers, upper, function(x, b) (b - x)/h)
        kernel <- .surface_normal_mass(a, b)
        if (metric == "thickness") {
            a <- outer(-centers, lower, function(x, b) (b - x)/h)
            b <- outer(-centers, upper, function(x, b) (b - x)/h)
            kernel <- kernel + .surface_normal_mass(a, b)
        }
        if (max(abs(rowSums(kernel) - 1)) > 1e-10) 
            stop("Kernel row mass is not one")
        .surface_kernel_cache[[tag]] <- list(klo = klo, khi = khi, kernel = kernel, breaks = cc$breaks)
    }
    ncenters <- khi - klo + 1L
    linear <- gg + (ix - klo) * length(rois)
    sums <- rowsum(ww, group = linear, reorder = FALSE)
    W <- matrix(0, length(rois), ncenters)
    W[as.integer(rownames(sums))] <- as.numeric(sums)
    P <- W %*% kernel
    if (any(P < -1e-12) || any(!is.finite(P)) || any(abs(rowSums(P)[valid] - 1) > 1e-10)) 
        stop("Fast KDE mass failed")
    P[P < 0] <- 0
    P[valid, ] <- P[valid, , drop = FALSE]/rowSums(P)[valid]
    P[!valid, ] <- NA_real_
    rownames(P) <- names(rois)
    P
}

surface_fast_masses_sparse <-
function (rois, metric, config, spacing_divisor = 32L) 
{
    cc <- config[[metric]]
    h <- cc$bandwidth
    delta <- h/spacing_divisor
    ix <- ww <- gg <- list()
    valid <- rep(FALSE, length(rois))
    nr <- length(rois)
    for (i in seq_along(rois)) {
        a <- rois[[i]]
        if (is.null(a) || !a$technical_pass || !isTRUE(a$metric_available[[metric]])) 
            next
        x <- a$x[[metric]]
        w <- a$area/sum(a$area)
        if (any(!is.finite(x)) || any(!is.finite(w))) 
            stop("Invalid raw data in sparse-center computation")
        lo <- floor(x/delta)
        fraction <- x/delta - lo
        ix[[i]] <- c(lo, lo + 1)
        ww[[i]] <- c(w * (1 - fraction), w * fraction)
        gg[[i]] <- rep(i, 2L * length(x))
        valid[i] <- TRUE
    }
    if (!any(valid)) 
        return(matrix(NA_real_, nr, length(cc$breaks) - 1L, dimnames = list(names(rois), NULL)))
    ix <- unlist(ix)
    ww <- unlist(ww)
    gg <- unlist(gg)
    centers <- sort(unique(ix))
    nc <- length(centers)
    sums <- rowsum(ww, group = gg + (match(ix, centers) - 1L) * nr, reorder = FALSE)
    W <- matrix(0, nr, nc)
    W[as.integer(rownames(sums))] <- as.numeric(sums)
    P <- matrix(0, nr, length(cc$breaks) - 1L)
    lower <- head(cc$breaks, -1L)
    upper <- tail(cc$breaks, -1L)
    for (start in seq.int(1L, nc, by = 256L)) {
        ii <- start:min(start + 255L, nc)
        x <- centers[ii] * delta
        K <- .surface_normal_mass(outer(x, lower, function(x, b) (b - x)/h), outer(x, upper, function(x, 
            b) (b - x)/h))
        if (metric == "thickness") 
            K <- K + .surface_normal_mass(outer(-x, lower, function(x, b) (b - x)/h), outer(-x, upper, 
                function(x, b) (b - x)/h))
        if (max(abs(rowSums(K) - 1)) > 1e-10) 
            stop("Sparse kernel mass is not one")
        P <- P + W[, ii, drop = FALSE] %*% K
    }
    if (any(P < -1e-12) || any(!is.finite(P)) || any(abs(rowSums(P)[valid] - 1) > 1e-10)) 
        stop("Sparse KDE mass failed")
    P[P < 0] <- 0
    P[valid, ] <- P[valid, , drop = FALSE]/rowSums(P)[valid]
    P[!valid, ] <- NA_real_
    rownames(P) <- names(rois)
    attr(P, "NumericMethod") <- "same_h32_linear_binning_sparse_occupied_centers"
    P
}

surface_jsd_bits <-
function (p, q) 
{
    if (!is.numeric(p) || !is.numeric(q) || length(p) != length(q) || !length(p)) 
        stop("Probability vectors must have equal positive length")
    if (any(!is.finite(c(p, q))) || any(c(p, q) < 0) || abs(sum(p) - 1) > 1e-10 || abs(sum(q) - 1) > 
        1e-10) 
        stop("Inputs must be probability masses summing to 1")
    if (all(p == q)) 
        return(0)
    log_m <- log2(p + q) - 1
    ip <- p > 0
    iq <- q > 0
    answer <- 0.5 * sum(p[ip] * (log2(p[ip]) - log_m[ip])) + 0.5 * sum(q[iq] * (log2(q[iq]) - log_m[iq]))
    if (!is.finite(answer) || answer < -1e-12 || answer > 1 + 1e-12) 
        stop("JSD outside theoretical bounds")
    max(0, min(1, answer))
}

surface_magic <-
function (con) 
{
    x <- readBin(con, "raw", n = 3L)
    if (length(x) != 3L) 
        stop("Truncated magic")
    sum(as.integer(x) * c(65536, 256, 1))
}

surface_matrix_from_mass <-
function (P) 
{
    n <- nrow(P)
    out <- matrix(NA_real_, n, n, dimnames = list(rownames(P), rownames(P)))
    valid <- which(rowSums(!is.finite(P)) == 0)
    if (any(abs(rowSums(P[valid, , drop = FALSE]) - 1) > 1e-10) || any(P[valid, , drop = FALSE] < 0)) 
        stop("Invalid probability matrix")
    entropy <- function(x) {
        a <- x
        a[a <= 0] <- 1
        -rowSums(x * log2(a))
    }
    H <- entropy(P[valid, , drop = FALSE])
    for (k in seq_along(valid)) {
        i <- valid[k]
        out[i, i] <- 0
        later <- which(valid > i)
        if (length(later)) {
            j <- valid[later]
            mix <- sweep(P[j, , drop = FALSE], 2, P[i, ], "+")/2
            d <- entropy(mix) - 0.5 * (H[k] + H[later])
            if (any(d < -1e-10 | d > 1 + 1e-10)) 
                stop("JSD outside bounds")
            d <- pmax(0, pmin(1, d))
            out[i, j] <- out[j, i] <- d
        }
    }
    out
}

surface_provenance <-
function (path) 
{
    readlog <- function(f) {
        p <- file.path(path, "scripts", f)
        if (file.exists(p)) 
            readLines(p, warn = FALSE)
        else character()
    }
    stamp <- readlog("build-stamp.txt")
    rec <- readlog("recon-all.log")
    fast <- readlog("recon-surf.done")
    version <- grep("^VERSION ", fast, value = TRUE)
    cmds <- grep("--curv-map|--thickness|mris_curvature|mris_thickness", rec, value = TRUE)
    fastlog <- readlog("recon-surf.log")
    cmdall <- unique(c(cmds, grep("--curv-map|--thickness|mris_curvature|mris_thickness", fastlog, value = TRUE)))
    whitecurv <- all(vapply(c("lh", "rh"), function(h) any(grepl(paste0(h, "[.]white.*", h, "[.]curv"), 
        cmdall)), logical(1)))
    done <- if (length(fastlog)) 
        file.exists(file.path(path, "scripts", "recon-surf.done"))
    else file.exists(file.path(path, "scripts", "recon-all.done"))
    stream <- if (grepl("[.]long[.]", basename(path)) || any(grepl("(^| )-long( |$)", c(rec, fastlog)))) 
        "longitudinal"
    else "cross_sectional_no_long_command_detected"
    list(Software = if (length(fastlog)) "FastSurfer_recon_surf_with_FreeSurfer" else "FreeSurfer", FastSurferVersion = if (length(version)) sub("^VERSION ", 
        "", version[1]) else NA_character_, FreeSurferBuild = paste(stamp, collapse = ";"), Stream = stream, 
        CompletionMarker = done, WhiteCurvatureCommandVerified = whitecurv, Commands = cmdall, SourceLog = tail(c(rec, 
            fastlog), 40L))
}

surface_read_annot <-
function (path) 
{
    con <- file(path, "rb")
    on.exit(close(con))
    nv <- surface_read_int(con)
    if (nv < 1 || nv > 1e+07) 
        stop("Invalid annotation vertex count")
    pairs <- matrix(surface_read_int(con, 2L * nv), ncol = 2L, byrow = TRUE)
    if (anyDuplicated(pairs[, 1]) || !identical(sort(pairs[, 1]), 0:(nv - 1L))) 
        stop("Annotation vertex indices not a complete unique set")
    if (surface_read_int(con) == 0) 
        stop("Annotation has no embedded color table")
    n <- surface_read_int(con)
    rows <- list()
    if (n > 0) {
        origin <- surface_read_string(con)
        for (i in seq_len(n)) {
            nm <- surface_read_string(con)
            rgba <- surface_read_int(con, 4L)
            rows[[i]] <- data.frame(Index = i - 1L, Name = nm, R = rgba[1], G = rgba[2], B = rgba[3], 
                A = rgba[4])
        }
    }
    else {
        if (n != -2L) 
            stop("Unsupported annotation color-table version")
        maxidx <- surface_read_int(con)
        origin <- surface_read_string(con)
        nread <- surface_read_int(con)
        if (nread < 1 || nread > maxidx || maxidx > 10000) 
            stop("Invalid color-table dimensions")
        for (i in seq_len(nread)) {
            idx <- surface_read_int(con)
            nm <- surface_read_string(con)
            rgba <- surface_read_int(con, 4L)
            rows[[i]] <- data.frame(Index = idx, Name = nm, R = rgba[1], G = rgba[2], B = rgba[3], A = rgba[4])
        }
    }
    ct <- do.call(rbind, rows)
    if (anyDuplicated(ct$Index) || anyDuplicated(ct$Name) || any(as.matrix(ct[, c("R", "G", "B", "A")]) < 
        0) || any(as.matrix(ct[, c("R", "G", "B", "A")]) > 255)) 
        stop("Invalid color table")
    ct$Code <- ct$R + 256 * ct$G + 65536 * ct$B
    if (anyDuplicated(ct$Code)) 
        stop("Ambiguous annotation colors")
    codes <- pairs[order(pairs[, 1]), 2]
    labels <- ct$Name[match(codes, ct$Code)]
    if (any(is.na(labels) & codes != 0L)) 
        stop("Unrecognized nonzero annotation color")
    labels[is.na(labels)] <- "unknown"
    list(labels = labels, codes = codes, color_table = ct, origin = origin)
}

surface_read_cortex <-
function (path, nv) 
{
    h <- readLines(path, n = 2L, warn = FALSE)
    if (length(h) != 2L) 
        stop("Truncated cortex label")
    n <- as.integer(h[2])
    d <- read.table(path, skip = 2L, header = FALSE, comment.char = "", colClasses = "numeric")
    if (nrow(d) != n || ncol(d) != 5L || anyDuplicated(d[, 1]) || any(d[, 1] < 0 | d[, 1] >= nv | d[, 
        1] != floor(d[, 1]))) 
        stop("Invalid cortex label indices")
    cortex <- rep(FALSE, nv)
    cortex[d[, 1] + 1L] <- TRUE
    cortex
}

surface_read_float <-
function (con, n) 
{
    x <- readBin(con, numeric(), n = n, size = 4L, endian = "big")
    if (length(x) != n) 
        stop("Truncated float data")
    x
}

surface_read_int <-
function (con, n = 1L) 
{
    x <- readBin(con, integer(), n = n, size = 4L, endian = "big")
    if (length(x) != n) 
        stop("Truncated integer data")
    x
}

surface_read_mesh <-
function (path) 
{
    con <- file(path, "rb")
    on.exit(close(con))
    if (surface_magic(con) != 16777214) 
        stop("Only triangular FreeSurfer meshes supported")
    header <- readLines(con, n = 2L, warn = FALSE)
    h <- surface_read_int(con, 2L)
    if (h[1] < 3 || h[1] > 1e+07 || h[2] < 1 || h[2] > 2e+07) 
        stop("Invalid mesh header")
    vertices <- matrix(surface_read_float(con, h[1] * 3L), ncol = 3L, byrow = TRUE)
    faces <- matrix(surface_read_int(con, h[2] * 3L), ncol = 3L, byrow = TRUE)
    list(vertices = vertices, faces = faces, header = header)
}

surface_read_morph <-
function (path) 
{
    con <- file(path, "rb")
    on.exit(close(con))
    if (surface_magic(con) != 16777215) 
        stop("Only modern one-value-per-vertex morphometry supported")
    h <- surface_read_int(con, 3L)
    if (h[1] < 1 || h[1] > 1e+07 || h[2] < 1 || h[3] != 1L) 
        stop("Invalid morphometry header")
    x <- surface_read_float(con, h[1])
    if (length(readBin(con, "raw", n = 1L))) 
        stop("Unexpected trailing morphometry data")
    attr(x, "faces") <- h[2]
    x
}

surface_read_string <-
function (con) 
{
    n <- surface_read_int(con)
    if (n < 1 || n > 1e+06) 
        stop("Invalid annotation string length")
    bytes <- readBin(con, "raw", n = n)
    if (length(bytes) != n) 
        stop("Truncated annotation string")
    if (tail(bytes, 1) == as.raw(0)) 
        bytes <- head(bytes, -1L)
    rawToChar(bytes)
}

surface_roi_names <-
function () 
{
    r <- c("caudalanteriorcingulate", "caudalmiddlefrontal", "cuneus", "entorhinal", "fusiform", "inferiorparietal", 
        "inferiortemporal", "isthmuscingulate", "lateraloccipital", "lateralorbitofrontal", "lingual", 
        "medialorbitofrontal", "middletemporal", "parahippocampal", "paracentral", "parsopercularis", 
        "parsorbitalis", "parstriangularis", "pericalcarine", "postcentral", "posteriorcingulate", "precentral", 
        "precuneus", "rostralanteriorcingulate", "rostralmiddlefrontal", "superiorfrontal", "superiorparietal", 
        "superiortemporal", "supramarginal", "transversetemporal", "insula")
    c(paste0("lh_", r), paste0("rh_", r))
}

surface_vertex_areas <-
function (vertices, faces_zero_based) 
{
    stopifnot(is.matrix(vertices), ncol(vertices) == 3L, all(is.finite(vertices)), is.matrix(faces_zero_based), 
        ncol(faces_zero_based) == 3L, all(is.finite(faces_zero_based)), all(faces_zero_based == floor(faces_zero_based)), 
        all(faces_zero_based >= 0), all(faces_zero_based < nrow(vertices)))
    if (!nrow(faces_zero_based)) 
        stop("No faces")
    faces <- faces_zero_based + 1L
    a <- vertices[faces[, 2], , drop = FALSE] - vertices[faces[, 1], , drop = FALSE]
    b <- vertices[faces[, 3], , drop = FALSE] - vertices[faces[, 1], , drop = FALSE]
    cross <- cbind(a[, 2] * b[, 3] - a[, 3] * b[, 2], a[, 3] * b[, 1] - a[, 1] * b[, 3], a[, 1] * b[, 
        2] - a[, 2] * b[, 1])
    fa <- sqrt(rowSums(cross^2))/2
    if (any(!is.finite(fa))) 
        stop("Nonfinite triangle area")
    out <- numeric(nrow(vertices))
    sums <- rowsum(rep(fa/3, times = 3L), group = as.vector(faces), reorder = FALSE)
    out[as.integer(rownames(sums))] <- as.numeric(sums)
    attr(out, "total_area") <- sum(fa)
    attr(out, "zero_area_faces") <- sum(fa == 0)
    out
}

surface_weighted_stats <-
function (x, w) 
{
    w <- w/sum(w)
    mu <- sum(w * x)
    s <- sqrt(sum(w * (x - mu)^2))
    o <- order(x)
    cum <- cumsum(w[o])
    q <- vapply(c(0.01, 0.05, 0.25, 0.5, 0.75, 0.95, 0.99), function(p) x[o][which(cum >= p)[1]], numeric(1))
    c(N = length(x), EffectiveN = 1/sum(w^2), Mean = mu, SD = s, Min = min(x), Max = max(x), P01 = q[1], 
        P05 = q[2], P25 = q[3], Median = q[4], P75 = q[5], P95 = q[6], P99 = q[7])
}

.numeric_matrix <-
function (x, name = "x") 
{
    x <- as.matrix(x)
    if (!is.numeric(x) || !length(x) || any(!is.finite(x))) 
        stop(name, " must be a nonempty finite numeric matrix")
    storage.mode(x) <- "double"
    x
}

