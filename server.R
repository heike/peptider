library(shiny)
library(peptider)
library(quantreg)
library(plyr)
library(dplyr)
library(scales)
library(ggplot2)
library(grid)
library(reshape2)

## Global Vars
scheme.list <- c("20/20" = "2020", "NNK/S" = "NNK", "NNB" = "NNB", "NNN" = "NNN", "20/20 (-C)" = "2020C", "NNK/S (-C)" = "NNKC", "NNB (-C)" = "NNBC", "NNN (-C)" = "NNNC")
k <- 1:20
N <- as.vector(sapply(10^seq(1, 25, by = 1), `*`, seq(1.0, 9.9, by = 0.1)))
degree <- 1:2
curscheme <- "2020"

## Global Data
load("data/library_initial.RData")
load("data/peptide_initial.RData")
load("data/neighborhood_initial.RData")

## Functions
fancy_scientific <- function(l) { 
    # turn in to character string in scientific notation 
    l <- format(l, scientific = TRUE) 
    # quote the part before the exponent to keep all the digits 
    l <- gsub("^(.*)e", "'\\1'e", l) 
    # turn the 'e+' into plotmath round 
    l <- gsub("e", "%*%10^", l) 
    # remove +
    l <- gsub("\\+", "", l)
    # return this as an expression 
    parse(text=l) 
}

specify_decimal <- function(x, k) format(round(x, k), nsmall=k, scientific = FALSE)

## Server Definition
shinyServer(function(input, output, session) {
    
    custom_aacids <- list()
    custom_cs <- list()
    
    custom_scheme <- reactive({
        if (input$custom_done == 0) return()

        isolate({
            aacid <- sapply(c(1:(length(custom_aacids) + 1), 0), function(i) {
                input[[paste("custom_aacid", i, sep = "")]]
            })
            aacid_clean <- aacid[aacid != ""]
            c <- sapply(c(1:(length(custom_aacids) + 1), 0), function(i) {
                input[[paste("custom_c", i, sep = "")]]
            })
            c_clean <- c[aacid != ""]
            class_clean <- c(LETTERS[(length(c_clean) - 1):1], "Z")
            
            my.df <- data.frame(class = class_clean, aacid = aacid_clean, c = c_clean)
            my.df$class <- factor(my.df$class, levels = c(LETTERS[1:(length(c_clean) - 1)], "Z"))
            
            return(arrange(my.df, class))
        })
    })
    
    output$custom_aacid <- renderUI({
        if (input$custom_new == 0) return()
        
        custom_aacids[[length(custom_aacids) + 1]] <<- textInput(paste("custom_aacid", length(custom_aacids) + 2, sep = ""), "", value="")
        
        custom_aacids
    })
    
    output$custom_c <- renderUI({
        if (input$custom_new == 0) return()
        
        custom_cs[[length(custom_cs) + 1]] <<- numericInput(paste("custom_c", length(custom_cs) + 2, sep = ""), "", value = length(custom_cs) + 2)
        
        custom_cs
    })

    ## Data    
    custom.data <- reactive({
        if (!is.null(custom_scheme())) {
            withProgress(message = 'Calculating, please wait...', detail = 'This may take up to a minute...', value = 0, {
                val <- peptider:::generateCustom(scheme_def = custom_scheme(), savefile = FALSE)
                
                return(val)
            })
        } else {
            return(NULL)
        }
    })
    
    library.data <- reactive({
        withProgress(message = 'Calculating, please wait...', detail = 'This may take up to a minute...', value = 0, {
            
            lib.data <- library.initial
            if (my.scheme() == "Custom") lib.data <- rbind(lib.data, custom.data()$custom.lib)
            
            return(filter(lib.data, scheme %in% c(my.scheme(), input$schemes), k %in% c(input$length, input$lengths)))
        })
    })
    
    library.mine <- reactive({
        return(filter(library.data(), scheme == my.scheme(), k == input$length, N == N.mine()))
    })
    
    peptide.data <- reactive({
        withProgress(message = 'Calculating, please wait...', detail = 'This may take up to a minute...', value = 0, {
            pep.data <- peptide.initial
            if (my.scheme() == "Custom") pep.data <- rbind(pep.data, custom.data()$custom.probs)
            
            return(filter(pep.data, scheme %in% c(my.scheme(), input$schemes), k %in% c(input$length, input$lengths)))
        })
    })
    
    peptide.mine <- reactive({
        return(filter(peptide.data(), scheme == my.scheme(), k == input$length))
    })
    
    my.peptide <- reactive({
        return(toupper(input$peptide))
    })
    
    my.scheme <- reactive({
        schm <- input$scheme
        if (!is.null(custom_scheme()) & schm == "Custom") {
            curscheme <<- "Custom"
        } else if (schm != "Custom") {
            curscheme <<- schm
        }
        
        return(curscheme)
    })
    
    my.schemedef <- reactive({
        schmdef <- NULL
        if (my.scheme() == "Custom") {
            schmdef <- custom.data()$scheme_def
        } else {
            schmdef <- scheme(my.scheme())
        }
        
        return(schmdef)
    })
    
    N.mine <- reactive({
        return(input$size_x * 10^input$size_y)
    })
    
    inclusion_size <- reactive({
        selectedRange <- 10^(seq(7, 11))
        nValues <- if(N.mine() %in% selectedRange) selectedRange else c(selectedRange, N.mine())
        
        inclusion_size <- ldply(nValues, function(n) {
            cbind(peptide.mine(), N = n, inclusion = with(peptide.mine(), 1 - exp(-n * probs / di)))
        })
        inclusion_size$N <- factor(inclusion_size$N)
        
        return(inclusion_size)
    })
    
    neighborhood_cases <- reactive({
        nei.data <- neighborhood.initial
        if (my.scheme() == "Custom") nei.data <- rbind(nei.data, custom.data()$custom.nei)
        
        return(filter(nei.data, scheme %in% c(my.scheme(), input$schemes), k == input$length))
    })
    
    neighborhood_cases_mine <- reactive({
        return(filter(neighborhood_cases(), N == N.mine(), k == input$length, scheme == my.scheme()))
    })
    
    neighborhood_peptide <- reactive({
        withProgress(message = 'Calculating, please wait...', detail = 'This may take up to a minute...', value = 0, {
            nei1 <- getNeighbors(my.peptide())
            nei2 <- unique(unlist(getNeighbors(nei1)))
            
            nei.dat <- ldply(unique(c(input$schemes, my.scheme())), function(schm) {
                myschm <- schm
                if (myschm == "Custom") myschm <- my.schemedef()
                data.frame(scheme = schm, N = N,
                           "1" = ppeptide(x = nei1, libscheme = myschm, N = N),
                           "2" = ppeptide(x = nei2, libscheme = myschm, N = N))
            })
            
            nei.melt <- melt(nei.dat, id.vars = c("scheme", "N"), variable.name = "degree")
            nei.melt$degree <- gsub("X", "", nei.melt$degree)
            nei.melt$scheme <- factor(nei.melt$scheme, levels = scheme.list)
        })
        
        return(nei.melt)
    })
    
    neighborhood_mine <- reactive({
        return(filter(neighborhood_peptide(), scheme == my.scheme(), N == N.mine()))
    })
    
    getEncodingClass <- function(scheme, peptide) {
        lib <- peptider:::libscheme(my.schemedef(), 1)$info$scheme
        my.classes <- lib$class[sapply(strsplit(peptide, "")[[1]], grep, x = lib$aacid)]
        
        return(paste(my.classes, collapse = ""))
    }
    
    getReducedEncodingClass <- function(scheme, peptide) {
        lib <- peptider:::libscheme(my.schemedef(), 1)$info$scheme
        my.classes <- lib$class[sapply(strsplit(peptide, "")[[1]], grep, x = lib$aacid)]
        my.classes <- factor(my.classes, levels = unique(lib$class[lib$class != "Z"]))
        
        return(paste(as.numeric(table(my.classes)), collapse = ","))
    }
    
    ## Summary
    
    output$properties_table <- renderDataTable({        
        libprops <- c(specify_decimal(c(library.mine()$coverage, library.mine()$efficiency, library.mine()$diversity), input$digits), format(peptider:::diversity(input$length, my.schemedef(), N.mine()), scientific = TRUE, digits = input$digits + 1))
        
        my.df <- data.frame(Scheme = my.scheme(), Size = N.mine(), "Pept. Length" = input$length, Coverage = libprops[1], Efficiency = libprops[2], "Pept. Diversity" = libprops[4], "Func. Diversity" = libprops[3], check.names = FALSE)
        
        return(my.df)
    }, options = list(paging = FALSE, searching = FALSE))
    
    output$pep_properties_table <- renderDataTable({
        withProgress(message = 'Calculating, please wait...', detail = 'This may take up to a minute...', value = 0, {
                
            inclusion.mine <- filter(inclusion_size(), N == N.mine())
            inclusion.pep <- filter(inclusion.mine, class == getReducedEncodingClass(my.scheme(), my.peptide()))
            
            incl.probs <- specify_decimal(c(inclusion.pep$inclusion, min(inclusion.mine$inclusion), max(inclusion.mine$inclusion)), input$digits)
            
            nei.mine <- neighborhood_mine()
            nei.cases <- neighborhood_cases_mine()
            nei.my.probs <- as.numeric(nei.mine$value)
            
            degree1.probs <- specify_decimal(c(nei.my.probs[1], nei.cases$worst[nei.cases$degree == 1], nei.cases$best[nei.cases$degree == 1]), input$digits)
            degree2.probs <- specify_decimal(c(nei.my.probs[2], nei.cases$worst[nei.cases$degree == 2], nei.cases$best[nei.cases$degree == 2]), input$digits)
            
            cdns <- c(codons(my.peptide(), my.schemedef()), min(my.schemedef()$c)^input$length, max(my.schemedef()$c)^input$length)
            
            if (nchar(my.peptide()) == input$length) my.df <- data.frame(Peptide = c(my.peptide(), "Worst Case", "Best Case"), "DNA Encodings" = cdns, "Incl. Probability" = incl.probs, "Neighborhood D1" = degree1.probs, "Neighborhood D2" = degree2.probs, check.names = FALSE)
            else my.df <- data.frame(Peptide = c("Worst Case", "Best Case"), "DNA Encodings" = cdns[2:3], "Incl. Probability" = incl.probs[2:3], "Neighborhood D1" = degree1.probs[2:3], "Neighborhood D2" = degree2.probs[2:3], check.names = FALSE)
                
            return(my.df)
        })
    }, options = list(paging = FALSE, searching = FALSE))
    
    output$scheme_table <- renderDataTable({
        my.df <- peptider:::libscheme(my.schemedef(), 1)$info$scheme
        names(my.df) <- c("aa Class", "aa (Amino Acids)", "Codons", "aa Count")

        return(my.df)
    }, options = list(paging = FALSE, searching = FALSE))
    
    output$sample_table <- renderDataTable({
        withProgress(message = 'Calculating, please wait...', detail = 'This may take up to a minute...', value = 0, {
            lib <- peptider:::libscheme(my.schemedef(), 1)$info$scheme
            possible.aa <- unlist(strsplit(as.character(lib$aacid[lib$class != "Z"]), ""))
            
            pep.table <- ldply(1:10, function(y) {
                c("Peptide" = paste(sample(possible.aa, size = input$length), collapse = ""))
            })
            
            pep.encs <- sapply(pep.table$Peptide, function(y){getReducedEncodingClass(my.scheme(), y)})
            
            pep.table[,"Peptide Class"] <- sapply(pep.table$Peptide, function(y){getEncodingClass(my.scheme(), y)})
            pep.table$`DNA Encodings` <- codons(pep.table$Peptide, my.schemedef())
            pep.table[,"Incl. Probability"] <- sapply(pep.encs, function(y){inclusion_size()$inclusion[inclusion_size()$class == y & inclusion_size()$N == N.mine()]})
                
            pep.table[,"Incl. Probability"] <- specify_decimal(pep.table[,"Incl. Probability"], input$digits)
            
            return(pep.table)
        })
    }, options = list(paging = FALSE, searching = FALSE))
    
    ## Inclusion
    output$inclusion_text <- renderText({
        inclusion.mine <- filter(inclusion_size(), N == N.mine())
        return(paste("The probability that a peptide is part of your library ranges from", specify_decimal(min(inclusion.mine$inclusion), input$digits), "to", specify_decimal(max(inclusion.mine$inclusion), input$digits)))
    })
    
    output$inclusion_text_pep <- renderText({
        inclusion.mine <- filter(inclusion_size(), N == N.mine(), class == getReducedEncodingClass(my.scheme(), my.peptide()))
        return(paste("The probability that", my.peptide(), "is part of your library is", specify_decimal(inclusion.mine$inclusion, input$digits)))
    })
    
    output$inclusion_caption <- renderText({
        return("The probability that a peptide is included in the defined library depends on the encoding scheme and its sequence. Here, the highest and lowest probabilities possible for peptide sequences of the library are determined. All peptide inclusion probabilities of peptides from the library are found within this range. This range also defines the chances that the \"best\" possible peptide is included in the library (\"best\", in this case, is defined as the peptides that meets the experimental goals best). For individual peptides or encoding schemes with a one codon/amino acid ratio specific probabilities can be given.")
    })
    
    output$inclusion_plot <- renderPlot({
        withProgress(message = 'Calculating, please wait...', detail = 'This may take up to a minute...', value = 0, {
            selectedRange <- 10^(seq(7, 11))
            nValues <- if(N.mine() %in% selectedRange) selectedRange else c(selectedRange, N.mine())
            inclusion.mine <- filter(inclusion_size(), N == N.mine(), class == getReducedEncodingClass(my.scheme(), my.peptide()))
            
            print(
                qplot(N, inclusion, geom = "boxplot", weight = di*choices, data = inclusion_size(), fill = N, group = N, alpha = I(0.8)) +
                    scale_x_discrete(breaks = factor(nValues), labels = c(fancy_scientific(1 * 10^(seq(7, 11))), fancy_scientific(N.mine()))[1:length(nValues)]) + 
                    scale_fill_brewer(palette="Set2", guide="none") + 
                    ylim(0, 1) +
                    ylab("Probability of Inclusion") +
                    xlab("Library Size (N)") +
                    coord_flip() +
                    theme_bw() +
                    theme(legend.position = "bottom") +
                    if (input$length == nchar(my.peptide())) {
                        list(geom_point(data = inclusion.mine, aes(x = N, y = inclusion, shape = k), col = I("black"), size = 6, inherit.aes = FALSE),
                        scale_shape_manual("", values = 1, labels = my.peptide()))
                    }
            )
        })
    })
    
    output$inclusion_caption_schemes <- renderText({
        return(paste("Overview of the influence of different encoding schemes on the probability of peptide enclosure within a library from", N.mine(), "DNA sequences and a peptide length of", input$length, "amino acids. The probability that", my.peptide(), "is in your library is displayed on the plot as a circle."))
    })
    
    output$inclusion_plot_schemes <- renderPlot({
        withProgress(message = 'Calculating, please wait...', detail = 'This may take up to a minute...', value = 0, {
            peptide.sub <- filter(peptide.data(), k == input$length)
            
            inclusion_schemes <- cbind(peptide.sub, inclusion = with(peptide.sub, 1 - exp(-N.mine() * probs / di)))
            inclusion.mine <- filter(inclusion_size(), N == N.mine(), class == getReducedEncodingClass(my.scheme(), my.peptide()))
            
            print(
                ggplot(inclusion_schemes, aes(scheme, inclusion)) + 
                    geom_boxplot(aes(weight = di * choices, fill = scheme), alpha = 0.8) +
                    scale_fill_brewer(palette="Set1", guide="none") +
                    ylim(0, 1) +
                    ylab("Probability of Inclusion") +
                    xlab("Library Scheme") +
                    coord_flip() +
                    theme_bw() +
                    theme(legend.position = "bottom") +
                    if (input$length == nchar(my.peptide())) {
                        list(geom_point(data = inclusion.mine, aes(x = scheme, y = inclusion, shape = k), col = I("black"), size = 6, inherit.aes = FALSE),
                             scale_shape_manual("", values = 1, labels = my.peptide()))
                    }
            )
        })
    })
    
    ## Neighborhoods
    output$neighborhood_text_degree1 <- renderText({        
        deg1 <- specify_decimal(as.numeric(neighborhood_cases_mine()[neighborhood_cases_mine()$degree == 1, c("worst", "best")]), input$digits)
        
        return(paste("The probability that a top performing peptide (degree 1) is encoded by your library ranges from", min(deg1), "to", max(deg1)))
    })
    
    output$neighborhood_text_degree2 <- renderText({
        deg2 <- specify_decimal(as.numeric(neighborhood_cases_mine()[neighborhood_cases_mine()$degree == 2, c("worst", "best")]), input$digits)
        
        return(paste("The probability that a top performing peptide (degree 2) is encoded by your library ranges from", min(deg2), "to", max(deg2)))
    })
    
    output$neighborhood_plot <- renderPlot({
        withProgress(message = 'Calculating, please wait...', detail = 'This may take up to a minute...', value = 0, {
            nei.melt <- melt(neighborhood_cases(), id.var = c("scheme", "N", "k", "degree"), variable.name="Scenario")
            nei.melt$scheme <- factor(nei.melt$scheme, levels = scheme.list)
            
            print(
                qplot(N, value, data = nei.melt, geom = "line", size = I(1.5), linetype = Scenario, colour = scheme) +
                    facet_grid(degree~., labeller=label_both) +
                    scale_x_log10(breaks=c(10^(1:25)),labels=expression(10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10, 10^11, 10^12, 10^13, 10^14, 10^15, 10^16, 10^17, 10^18, 10^19, 10^20, 10^21, 10^22, 10^23, 10^24, 10^25)) +
                    scale_colour_brewer("Library\nScheme", palette="Set1") +
                    geom_vline(xintercept = library.mine()$N) +
                    geom_point(inherit.aes = FALSE, data = neighborhood_cases_mine(), aes(x = N, y = worst, shape = "2"), size = 6, colour = "black") +
                    geom_point(inherit.aes = FALSE, data = neighborhood_cases_mine(), aes(x = N, y = best, shape = "2"), size = 6, colour = "black") +
                    scale_shape_manual("", values = c(2, 2), labels = "Your Library") +
                    theme_bw()
            )
        })
    })
    
    output$neighborhood_caption <- renderText({
        return("Overview of the probability that at least one of the sequences belonging to the neighborhoods of a peptide sequence of degree one and two, respectively (varying from the reference sequence in up to one or up to two amino acid positions) is included in the library. This range also defines the probability that at least one \"top\" peptides is included in the library. A \"top\" peptide is defined as a peptide belonging to the degree 1 or 2 neighborhood of the \"best\" possible peptide (see \"Inclusion\" tab). Best and worst case probabilities depend on the number of encodings for a sequence and the probabilities of exchangeability of amino acids.")
    })
    
    output$neighborhood_text_mine_degree1 <- renderText({
        nei.mine <- neighborhood_mine()
        
        return(paste("The probability that a peptide in the degree 1 neighborhood of", my.peptide(), "is part of your library is", specify_decimal(as.numeric(nei.mine[nei.mine$degree == 1, "value"]), input$digits)))
    })
    
    output$neighborhood_text_mine_degree2 <- renderText({       
        nei.mine <- neighborhood_mine()
        
        return(paste("The probability that a peptide in the degree 2 neighborhood of", my.peptide(), "is part of your library is", specify_decimal(as.numeric(nei.mine[nei.mine$degree == 2, "value"]), input$digits)))
    })
    
    output$neighborhood_plot_mine <- renderPlot({
        withProgress(message = 'Calculating, please wait...', detail = 'This may take up to a minute...', value = 0, {
            
            nei.mine <- neighborhood_peptide()
                    
            print(
                qplot(N, value, data = nei.mine, colour = scheme, geom = "line", size = I(1.5)) + 
                    facet_grid(degree~., labeller=label_both) +
                    scale_x_log10(breaks=c(10^(1:25)),labels=expression(10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10, 10^11, 10^12, 10^13, 10^14, 10^15, 10^16, 10^17, 10^18, 10^19, 10^20, 10^21, 10^22, 10^23, 10^24, 10^25)) +
                    ylab("Probability of at least one peptide from the neighborhood") + 
                    scale_colour_brewer("Library\nScheme", palette="Set1") + 
                    xlab("Library Size (N)") +
                    geom_vline(xintercept = library.mine()$N) +
                    geom_point(inherit.aes = FALSE, data = neighborhood_mine(), aes(x = N, y = value, shape = "2"), size = 6, colour = "black") +
                    scale_shape_manual("", values = c(2, 2), labels = "Your Library") +
                    theme_bw()    
            )
        })
    })
    
    output$neighborhood_caption_mine <- renderText({
        return(paste("Overview of the probability that at least one of the sequences belonging to the neighborhoods of", my.peptide(), "of degree one and two, respectively, is part of your library."))
    })
    
    ## Coverage
    output$coverage_text <- renderText({
        return(paste("The expected coverage of your library is", specify_decimal(library.mine()$coverage, input$digits), "+/-", format(sqrt(peptider:::coverage(input$length, my.schemedef(), N.mine(), variance = TRUE)), digits = input$digits + 1)))
    })
    
    output$coverage_caption <- renderText({
        return("Overview of expected coverage for selected libraries. The expected coverage of your library is shown on the plot as a triangle.")
    })
    
    output$coverage_plot <- renderPlot({
        withProgress(message = 'Calculating, please wait...', detail = 'This may take up to a minute...', value = 0, {
            print(
                qplot(N, coverage, data = library.data(), geom = "line", colour = scheme, linetype = k, size = I(1.2)) +
                geom_vline(xintercept = library.mine()$N) +
                xlab("Library Size (N)") + ylab("Coverage") +
                scale_x_log10(breaks=c(10^(1:25)),labels=expression(10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10, 10^11, 10^12, 10^13, 10^14, 10^15, 10^16, 10^17, 10^18, 10^19, 10^20, 10^21, 10^22, 10^23, 10^24, 10^25)) +
                scale_linetype("Sequence\nLength (k)") +  
                scale_colour_brewer("Library\nScheme", palette="Set1") + 
                theme_bw()  + 
                theme(legend.key.width = unit(1.5, "cm")) +
                geom_point(inherit.aes = FALSE, data = library.mine(), aes(x = N, y = coverage, shape = k), size = 6, colour = "black") +
                scale_shape_manual("", values = 2, labels = "Your Library")
            )  
        })
    })
    
    ## Efficiency
    output$efficiency_text <- renderText({
        return(paste("The relative efficiency of your library is", specify_decimal(library.mine()$efficiency, input$digits), "+/-", format(sqrt(peptider:::efficiency(input$length, my.schemedef(), N.mine(), variance = TRUE)), digits = input$digits + 1)))
    })
    
    output$efficiency_caption <- renderText({
        return("Overview of relative efficiency for selected libraries. The relative efficiency of your library is shown on the plot as a triangle. If library schemes are used that encode invalid peptides (e.g. containing stop codons) the maximum achievable relative efficiency drops below 1 (initial loss).")
    })
    
    output$efficiency_plot <- renderPlot({
        withProgress(message = 'Calculating, please wait...', detail = 'This may take up to a minute...', value = 0, {
            print(
                qplot(N, efficiency, data = library.data(), geom = "line", colour = scheme, linetype = k, size = I(1.2)) +
                    geom_vline(xintercept = library.mine()$N) +
                    xlab("Library Size (N)") + ylab("Relative Efficiency") +
                    ylim(c(0,1)) +
                    scale_x_log10(breaks=c(10^(1:25)),labels=expression(10^1, 10^2, 10^3, 10^4, 10^5, 10^6, 10^7, 10^8, 10^9, 10^10, 10^11, 10^12, 10^13, 10^14, 10^15, 10^16, 10^17, 10^18, 10^19, 10^20, 10^21, 10^22, 10^23, 10^24, 10^25)) +
                    scale_linetype("Sequence\nLength (k)") +  
                    scale_colour_brewer("Library\nScheme", palette="Set1") +
                    theme_bw()  + 
                    theme(legend.key.width = unit(1.5, "cm")) +
                    geom_point(inherit.aes = FALSE, data = library.mine(), aes(x = N, y = efficiency, shape = k), size = 6, colour = "black") +
                    scale_shape_manual("", values = 2, labels = "Your Library")
            )
        })
    })
    
    ## Diversity
    output$diversity_text <- renderText({
        return(paste("The functional diversity of your library is", specify_decimal(library.mine()$diversity, input$digits)))
    })
    
    output$pepdiversity_text <- renderText({
        return(paste("The peptide diversity of your library is", format(peptider:::diversity(input$length, my.schemedef(), N.mine()), scientific = TRUE, digits = input$digits + 1), "+/-", format(sqrt(peptider:::diversity(input$length, my.schemedef(), N.mine(), variance = TRUE)), scientific = TRUE, digits = input$digits + 1)))
    })
    
    output$pepdiversity_caption <- renderText({
        return("The peptide diversity is the number of distinct peptides encoded in the defined library")
    })
    
    output$diversity_caption <- renderText({
        return("Overview of functional diversity as defined by Makowski and colleges. The functional diversity does not describe the number of peptides found in a library, but can be used to describe the quality of the encoding scheme.")
    })
    
    output$diversity_table <- renderTable({
        diversity.table <- library.data() %>% group_by(scheme, k) %>% summarise(diversity = diversity[1])
        diversity.table$diversity <- specify_decimal(diversity.table$diversity, input$digits)
        
        diversity.dcast <- dcast(diversity.table, k ~ scheme, value.var = "diversity")
        names(diversity.dcast) <- c("Sequence Length (k)", names(scheme.list[scheme.list %in% c(input$schemes, my.scheme())]))
                
        return(diversity.dcast)
    }, include.rownames = FALSE)
    
})
