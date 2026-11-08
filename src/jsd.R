cortical_mri_stats <-
function (path) 
{
    a <- adni_extract_volumes(path, read_label_volumes = FALSE)
    roi <- surface_roi_names()
    list(cortex_volume_mm3 = unname(a$values["CortexVol_stats_mm3"]), roi = data.frame(ROI = roi, AREA_MM2 = unname(a$values[paste0("DKT62Area_stats_mm2__", 
        roi)]), VOLUME_MM3 = unname(a$values[paste0("DKT62Volume_stats_mm3__", roi)])), state = a$state)
}

edge_dictionary <-
function () 
{
    r <- surface_roi_names()
    ij <- which(upper.tri(matrix(FALSE, length(r), length(r))), arr.ind = TRUE)
    data.frame(EDGE_ID = seq_len(nrow(ij)), FEATURE = sprintf("JSD_%04d", seq_len(nrow(ij))), ROI1 = r[ij[, 
        1]], ROI2 = r[ij[, 2]], I = ij[, 1], J = ij[, 2])
}

extract_thickness_scan <-
function (path, config = thickness_config(), return_mass = FALSE) 
{
    a <- thickness_features(path, config, return_mass)
    s <- cortical_mri_stats(path)
    j <- match(config$roi_names, a$roi_stats$ROI)
    mu <- a$roi_stats$Mean[j]
    area <- s$roi$AREA_MM2[match(config$roi_names, s$roi$ROI)]
    avg <- if (all(is.finite(c(mu, area))) && sum(area) > 0) 
        global_mean_thickness(matrix(mu, nrow = 1), matrix(area, nrow = 1))
    else NA_real_
    a$cortical_stats <- s
    a$global_mri <- data.frame(CORTEX_CM3 = s$cortex_volume_mm3/1000, DKT62_AREA_WEIGHTED_THICKNESS_MM = avg)
    a
}

global_mean_thickness <-
function (means, areas) 
{
    means <- .numeric_matrix(means, "means")
    areas <- .numeric_matrix(areas, "areas")
    if (!identical(dim(means), dim(areas)) || ncol(means) != 62L || any(areas < 0) || any(rowSums(areas) <= 
        0) || !identical(colnames(means), colnames(areas))) 
        stop("Require matching DKT62 means and stats areas in the fixed ROI order")
    rowSums(means * areas)/rowSums(areas)
}

jsd_bits <-
function (p, q) 
surface_jsd_bits(p, q)

roi_dictionary <-
function () 
{
    data.frame(ROI_ID = seq_len(62), ROI = surface_roi_names(), HEMISPHERE = rep(c("lh", "rh"), each = 31), 
        stringsAsFactors = FALSE)
}

thickness_features <-
function (scan, config = thickness_config(), return_mass = FALSE) 
{
    raw <- if (is.character(scan) && length(scan) == 1L) 
        surface_extract_scan(scan)
    else scan
    if (!is.list(raw) || is.null(raw$ROI) || !identical(names(raw$ROI), config$roi_names)) 
        stop("Expected reconstructed scan directory or audited DKT62 raw cache")
    mass <- tryCatch(surface_fast_masses(raw$ROI, "thickness", config, config$linear_binning_divisor), 
        error = function(e) {
            if (!grepl("Extreme data range requires technical review", conditionMessage(e), fixed = TRUE)) 
                stop(e)
            surface_fast_masses_sparse(raw$ROI, "thickness", config, config$linear_binning_divisor)
        })
    mat <- surface_matrix_from_mass(mass)
    ed <- edge_dictionary()
    edges <- ed
    edges$JSD <- mat[cbind(ed$I, ed$J)]
    moments <- lapply(raw$ROI, function(a) {
        if (is.null(a) || !isTRUE(a$technical_pass) || !isTRUE(a$metric_available[["thickness"]])) 
            return(stats::setNames(rep(NA_real_, 13L), c("N", "EffectiveN", "Mean", "SD", "Min", "Max", 
                "P01", "P05", "P25", "Median", "P75", "P95", "P99")))
        weighted_roi_stats(a$x$thickness, a$area)
    })
    tab <- as.data.frame(do.call(rbind, moments))
    tab$ROI <- names(raw$ROI)
    result <- list(edges = edges, matrix = mat, roi_stats = tab, qc = raw$QC, provenance = raw$Provenance, 
        config = config, manual_anatomical_qc = "must be documented by the investigator")
    if (return_mass) 
        result$mass <- mass
    result
}

weighted_roi_stats <-
function (x, weights) 
{
    if (!is.numeric(x) || !length(x) || any(!is.finite(x)) || length(x) != length(weights) || any(!is.finite(weights)) || 
        any(weights < 0) || !is.finite(sum(weights)) || sum(weights) <= 0) 
        stop("Invalid measurements or weights")
    surface_weighted_stats(x, weights)
}


thickness_config <- function() dget(file.path(REPO_ROOT,"config","jsd.R"))
