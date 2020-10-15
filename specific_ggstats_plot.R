ggbetweenstats <- function(data,
                           x,
                           y,
                           plot.type = "box",
                           type = "parametric",
                           pairwise.comparisons = FALSE,
                           pairwise.annotation = "p.value",
                           pairwise.display = "significant",
                           p.adjust.method = "holm",
                           effsize.type = "unbiased",
                           partial = TRUE,
                           effsize.noncentral = TRUE,
                           bf.prior = 0.707,
                           bf.message = F,
                           results.subtitle = F,
                           xlab = NULL,
                           ylab = NULL,
                           caption = NULL,
                           title = NULL,
                           subtitle = NULL,
                           stat.title = NULL,
                           sample.size.label = F,
                           k = 2,
                           var.equal = FALSE,
                           conf.level = 0.95,
                           nboot = 100,
                           tr = 0.1,
                           mean.plotting = TRUE,
                           mean.ci = FALSE,
                           mean.point.args = list(size = 5, color = "darkred"), ### 15 before
                           mean.label.args = list(size = 5),
                           notch = FALSE,
                           notchwidth = 0.5,
                           linetype = "solid",
                           outlier.tagging = FALSE,
                           outlier.label = NULL,
                           outlier.coef = 1.5,
                           outlier.shape = 19,
                           outlier.color = "black",
                           outlier.label.args = list(size = 3),
                           outlier.point.args = list(),
                           point.args = list(
                             position = ggplot2::position_jitterdodge(dodge.width = 0.80),
                             alpha = 0.5,
                             size = 3,
                             stroke = 0,jitter.width = 0.6
                           ),
                           violin.args = list(width = 0.65, alpha = 0.2),
                           ggtheme = ggplot2::theme_bw(),
                           ggstatsplot.layer = TRUE,
                           package = "RColorBrewer",
                           palette = "Dark2",
                           direction = 1,
                           ggplot.component = NULL,
                           output = "plot",
                           messages = F,
                           ...) {
  
  # convert entered stats type to a standard notation
  type <- stats_type_switch(type)
  
  # no pairwise comparisons are available for Bayesian t-tests
  if (type == "bayes") pairwise.comparisons <- FALSE
  
  # ------------------------------ variable names ----------------------------
  
  # ensure the variables work quoted or unquoted
  x <- rlang::ensym(x)
  y <- rlang::ensym(y)
  outlier.label <- if (!rlang::quo_is_null(rlang::enquo(outlier.label))) {
    rlang::ensym(outlier.label)
  }
  
  # if `xlab` and `ylab` is not provided, use the variable `x` and `y` name
  if (is.null(xlab)) xlab <- rlang::as_name(x)
  if (is.null(ylab)) ylab <- rlang::as_name(y)
  
  # --------------------------------- data -----------------------------------
  
  # creating a dataframe
  data %<>%
    dplyr::select(.data = ., {{ x }}, {{ y }}, outlier.label = {{ outlier.label }}) %>%
    tidyr::drop_na(data = .) %>%
    dplyr::mutate(.data = ., {{ x }} := droplevels(as.factor({{ x }}))) %>%
    as_tibble(x = .)
  
  # if outlier.label column is not present, just use the values from `y` column
  if (rlang::quo_is_null(rlang::enquo(outlier.label))) {
    data %<>% dplyr::mutate(.data = ., outlier.label = {{ y }})
  }
  
  # add a logical column indicating whether a point is or is not an outlier
  data %<>%
    ipmisc::outlier_df(
      data = .,
      x = {{ x }},
      y = {{ y }},
      outlier.coef = outlier.coef,
      outlier.label = outlier.label
    )
  
  # figure out which test to run based on the number of levels of the
  # independent variables
  test <- ifelse(nlevels(data %>% dplyr::pull({{ x }}))[[1]] < 3, "t", "anova")
  
  # --------------------- subtitle/caption preparation ------------------------
  
  if (isTRUE(results.subtitle)) {
    # preparing the Bayes factor message
    if (type == "parametric" && isTRUE(bf.message)) {
      caption <-
        caption_function_switch(
          test = test,
          data = data,
          x = rlang::as_string(x),
          y = rlang::as_string(y),
          bf.prior = bf.prior,
          caption = caption,
          paired = FALSE,
          output = "caption",
          k = k
        )
    }
    
    # extracting the subtitle using the switch function
    subtitle <-
      subtitle_function_switch(
        # switch based on
        type = type,
        test = test,
        # arguments relevant for subtitle helper functions
        data = data,
        x = {{ x }},
        y = {{ y }},
        paired = FALSE,
        effsize.type = effsize.type,
        partial = partial,
        effsize.noncentral = effsize.noncentral,
        var.equal = var.equal,
        bf.prior = bf.prior,
        tr = tr,
        nboot = nboot,
        conf.level = conf.level,
        stat.title = stat.title,
        k = k,
        messages = messages
      )
  }
  
  # quit early if only subtitle is needed
  if (output %in% c("subtitle", "caption")) {
    return(switch(
      EXPR = output,
      "subtitle" = subtitle,
      "caption" = caption
    ))
  }
  
  # -------------------------- basic plot -----------------------------------
  
  # create the basic plot
  # add only the points which are *not* outliers
  plot <-
    ggplot2::ggplot(data = data, mapping = ggplot2::aes(x = {{ x }}, y = {{ y }})) +
    rlang::exec(
      .fn = ggplot2::geom_point,
      data = dplyr::filter(.data = data, !isanoutlier),
      na.rm = TRUE,
      ggplot2::aes(color = {{ x }}),
      !!!point.args
    )
  
  # if outliers are not being tagged, then add the points that were left out
  if (isFALSE(outlier.tagging)) {
    plot <- plot +
      rlang::exec(
        .fn = ggplot2::geom_point,
        data = dplyr::filter(.data = data, isanoutlier),
        na.rm = TRUE,
        ggplot2::aes(color = {{ x }}),
        !!!point.args
      )
  }
  
  # if outlier tagging is happening, decide how those points should be displayed
  if (isTRUE(outlier.tagging)) {
    if (plot.type == "violin") {
      plot <- plot +
        # add all outliers in
        ggplot2::geom_point(
          data = dplyr::filter(.data = data, isanoutlier),
          size = 3,
          stroke = 0,
          alpha = 0.7,
          na.rm = TRUE,
          color = outlier.color,
          shape = outlier.shape
        )
    }
  }
  
  # adding a boxplot
  if (plot.type %in% c("box", "boxviolin")) {
    if (isTRUE(outlier.tagging)) {
      plot <- plot +
        ggplot2::stat_boxplot(
          notch = notch,
          notchwidth = notchwidth,
          linetype = linetype,
          geom = "boxplot",
          width = 0.6,
          alpha = 0.2,
          fill = "white",
          outlier.shape = outlier.shape,
          outlier.size = 3,
          outlier.alpha = 0.7,
          outlier.color = outlier.color,
          coef = outlier.coef,
          na.rm = TRUE
        )
    } else {
      plot <- plot +
        ggplot2::geom_boxplot(
          notch = notch,
          notchwidth = notchwidth,
          linetype = linetype,
          width = 0.55,
          alpha = 0.2,
          fill=c("#FF0033","#3300FF"),
          outlier.shape = NA,
          position = ggplot2::position_dodge(width = NULL),
          na.rm = TRUE
        )
    }
  }
  
  # add violin geom
  if (plot.type %in% c("violin", "boxviolin")) {
    plot <- plot +
      rlang::exec(
        .fn = ggplot2::geom_violin,
        fill = "white",
        na.rm = TRUE,
        !!!violin.args
      )
  }
  
  # ---------------------------- outlier tagging -----------------------------
  
  # If `outlier.label` is not provided, outlier labels will just be values of
  # the `y` vector. If the outlier tag has been provided, just use the dataframe
  # already created.
  
  if (isTRUE(outlier.tagging)) {
    # applying the labels to tagged outliers with `ggrepel`
    plot <- plot +
      rlang::exec(
        .fn = ggrepel::geom_label_repel,
        data = dplyr::filter(.data = data, isanoutlier),
        mapping = ggplot2::aes(x = {{ x }}, y = {{ y }}, label = outlier.label),
        show.legend = FALSE,
        min.segment.length = 0,
        inherit.aes = FALSE,
        na.rm = TRUE,
        !!!outlier.label.args
      )
  }
  
  # ---------------- mean value tagging -------------------------------------
  
  # computing mean and confidence interval for mean using helper function
  # creating label column based on whether just mean is to be displayed or
  # mean plus its CI
  mean_dat <-
    mean_labeller(
      data = data,
      x = {{ x }},
      y = {{ y }},
      mean.ci = mean.ci,
      k = k
    )
  
  # add labels for mean values
  if (isTRUE(mean.plotting)) {
    plot <-
      mean_ggrepel(
        mean.data = mean_dat,
        x = {{ x }},
        y = {{ y }},
        plot = plot,
        mean.point.args = mean.point.args,
        mean.label.args = mean.label.args,
        inherit.aes = TRUE
      )
  }
  
  # ----------------- sample size labels --------------------------------------
  
  # adding sample size labels to the x axes
  if (isTRUE(sample.size.label)) {
    plot <- plot + ggplot2::scale_x_discrete(labels = c(unique(mean_dat$n_label)))
  }
  
  # ggsignif labels -----------------------------------------------------------
  
  if (isTRUE(pairwise.comparisons) && test == "anova") {
    # creating dataframe with pairwise comparison results
    df_pairwise <-
      pairwiseComparisons::pairwise_comparisons(
        data = data,
        x = {{ x }},
        y = {{ y }},
        type = type,
        tr = tr,
        paired = FALSE,
        var.equal = var.equal,
        p.adjust.method = p.adjust.method,
        k = k,
        messages = FALSE
      )
    
    # display the results if needed
    if (isTRUE(messages)) print(dplyr::select(df_pairwise, -label))
    
    # adding the layer for pairwise comparisons
    plot <-
      ggsignif_adder(
        plot = plot,
        df_pairwise = df_pairwise,
        data = data,
        x = {{ x }},
        y = {{ y }},
        pairwise.annotation = pairwise.annotation,
        pairwise.display = pairwise.display
      )
    
    # preparing the caption for pairwise comparisons test
    caption <- pairwise_caption(caption, unique(df_pairwise$test.details), p.adjust.method)
  }
  
  # ------------------------ annotations and themes -------------------------
  
  # specifying annotations and other aesthetic aspects for the plot
  plot <-
    aesthetic_addon(
      plot = plot,
      x = data %>% dplyr::pull({{ x }}),
      xlab = xlab,
      ylab = ylab,
      title = title,
      subtitle = subtitle,
      caption = caption,
      ggtheme = ggtheme,
      ggstatsplot.layer = ggstatsplot.layer,
      package = package,
      palette = palette,
      direction = direction,
      ggplot.component = ggplot.component
    )
  
  # --------------------- messages ------------------------------------------
  
  if (isTRUE(messages)) {
    # display normality test result as a message
    normality_message(
      x = data %>% dplyr::pull({{ y }}),
      lab = ylab,
      k = k
    )
    
    # display homogeneity of variance test as a message
    bartlett_message(
      data = data,
      x = {{ x }},
      y = {{ y }},
      lab = xlab,
      k = k
    )
  }
  
  # return the final plot
  return(plot)
}
stats_type_switch <- function(type) {
  dplyr::case_when(
    grepl("^p", type, TRUE) ~ "parametric",
    grepl("^n", type, TRUE) ~ "nonparametric",
    grepl("^r", type, TRUE) ~ "robust",
    grepl("^b", type, TRUE) ~ "bayes",
    TRUE ~ "parametric"
  )
}
caption_function_switch <- function(test, ...) {
  # choosing the appropriate test
  if (test == "t") {
    .f <- statsExpressions::bf_ttest
  } else {
    .f <- statsExpressions::bf_oneway_anova
  }
  
  # preparing the BF message for null
  rlang::exec(.fn = .f, ...)
}
subtitle_function_switch <- function(test, type, ...) {
  # figuring out type of test needed to run
  type <- stats_type_switch(type)
  
  # make a function character string
  .f_string <- paste("statsExpressions::expr_", test, "_", type, "(...)", sep = "")
  
  # evaluate it
  return(rlang::eval_bare(rlang::parse_expr(.f_string)))
}
mean_labeller <- function(data,
                          x,
                          y,
                          mean.ci = FALSE,
                          k = 3L,
                          ...) {
  
  # creating the dataframe
  data %<>%
    dplyr::select(.data = ., {{ x }}, {{ y }}) %>%
    tidyr::drop_na(data = .) %>%
    dplyr::mutate(.data = ., {{ x }} := droplevels(as.factor({{ x }}))) %>%
    as_tibble(x = .)
  
  # computing mean and confidence interval for mean
  mean_dat <-
    groupedstats::grouped_summary(
      data = data,
      grouping.vars = {{ x }},
      measures = {{ y }}
    ) %>% # introduce non-syntactic names to allow for `mean` pattern names
    dplyr::rename_at(
      .tbl = .,
      .vars = dplyr::vars(dplyr::matches("mean|^n$")),
      .funs = ~ paste(., "...summary", sep = "")
    ) %>%
    dplyr::mutate(.data = ., {{ y }} := `mean...summary`) %>%
    dplyr::select(.data = ., {{ x }}, {{ y }}, dplyr::contains("...")) %>%
    dplyr::mutate_at(
      .tbl = .,
      .vars = dplyr::vars(dplyr::matches("^mean\\.\\.\\.|^mean\\.conf")),
      .funs = ~ specify_decimal_p(x = ., k = k)
    ) %>%
    dplyr::group_nest(.tbl = ., {{ x }})
  
  # adding confidence intervals to the label for mean
  mean_dat %<>%
    dplyr::mutate(
      .data = .,
       label = data %>% {
         if (isTRUE(mean.ci)) {
           purrr::map(
             .x = .,
             .f = ~ paste(
                "list(~italic(widehat(mu))==",
                .$`mean...summary`,
               "",
                "CI[95*'%']",
                "*'['*",
                .$`mean.conf.low...summary`,
                ",",
                .$`mean.conf.high...summary`,
                "*']')",
               sep = ""
             )
           )
         } else {
           purrr::map(
             .x = .,
             .f = ~ paste("list(~italic(widehat(mu))==", .$`mean...summary`, ")", sep = " ")
           )
         }
       }
    )
  
  # adding sample size labels and arranging by original factor levels
   mean_dat %<>%
     tidyr::unnest(data = ., cols = c(data,label)) %>%
     dplyr::mutate(
       .data = .,
       n_label = paste0({{ x }}, "\n(n = ", `n...summary`, ")", sep = "")
     ) %>%
     dplyr::arrange(.data = ., {{ x }})

  # return the dataframe with mean information
  return(dplyr::select(mean_dat, -dplyr::contains("...")))
}
mean_ggrepel <- function(plot,
                         x,
                         y,
                         mean.data,
                         mean.point.args = list(size = 5, color = "darkred"),
                         mean.label.args = list(size = 3),
                         inherit.aes = TRUE,
                         ...) {
  # highlight the mean of each group
  plot <- plot +
    rlang::exec(
      .fn = ggplot2::stat_summary,
      mapping = ggplot2::aes(x = {{ x }}, y = {{ y }}),
      fun = mean,
      geom = "point",
      inherit.aes = inherit.aes,
      na.rm = TRUE,
      !!!mean.point.args
    )
  
  # attach the labels with means to the plot
  # plot +
  #   rlang::exec(
  #     .fn = ggrepel::geom_label_repel,
  #     data = mean.data,
  #     mapping = ggplot2::aes(x = {{ x }}, y = {{ y }}, label = label),
  #     show.legend = FALSE,
  #     min.segment.length = 0,
  #     inherit.aes = FALSE,
  #     parse = TRUE,
  #     na.rm = TRUE,
  #     !!!mean.label.args
  #   )
}
aesthetic_addon <- function(plot,
                            x,
                            xlab = NULL,
                            ylab = NULL,
                            title = NULL,
                            subtitle = NULL,
                            caption = NULL,
                            ggtheme = ggplot2::theme_bw(),
                            ggstatsplot.layer = TRUE,
                            package = "RColorBrewer",
                            palette = "Dark2",
                            direction = 1,
                            ggplot.component = NULL,
                            ...) {
  
  # if no. of factor levels is greater than the default palette color count
  palette_message(
    package = package,
    palette = palette,
    min_length = length(unique(levels(x)))[[1]]
  )
  
  # modifying the plot
  plot <- plot +
    ggplot2::labs(
      x = xlab,
      y = ylab,
      title = title,
      subtitle = subtitle,
      caption = caption,
      color = xlab
    ) +
    ggstatsplot::theme_ggstatsplot(
      ggtheme = ggtheme,
      ggstatsplot.layer = ggstatsplot.layer
    ) +
    ggplot2::theme(legend.position = "none") +
    paletteer::scale_color_paletteer_d(
      palette = paste0(package, "::", palette),
      direction = direction
    ) +
    paletteer::scale_fill_paletteer_d(
      palette = paste0(package, "::", palette),
      direction = direction
    )
  
  # ---------------- adding ggplot component ---------------------------------
  
  # return with any additional modification that needs to be made to the plot
  return(plot + ggplot.component)
}
palette_message <- function(package, palette, min_length) {
  # computing the number of colors in a given palette
  palette_df <-
    as_tibble(paletteer::palettes_d_names) %>%
    dplyr::filter(.data = ., package == !!package, palette == !!palette) %>%
    dplyr::select(.data = ., length)
  
  # if insufficient number of colors are available in a given palette
  if (palette_df$length[[1]] < min_length) {
    # message to display
    message(cat(
      ipmisc::red("Warning: "),
      ipmisc::blue("No. of factor levels is greater than default palette color count.\n"),
      ipmisc::blue("Try using another color `palette` (and/or `package`).\n")
    ),
    sep = ""
    )
  }
}