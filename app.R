# ============================================================
# Drug Sensitivity Explorer: Real GDSC2 Data
# Author: Gokul Selvaraj
# GitHub: GokulSelvaraj-Scientist
# Description: Interactive Shiny dashboard for exploring
#              drug sensitivity across cancer types using
#              real data from the Genomics of Drug Sensitivity
#              in Cancer (GDSC2) database
# Data: GDSC2 release 8.5 (Sanger Institute)
#       286 drugs, 969 cancer cell lines, 242,036 measurements
# ============================================================

# install.packages(c("shiny","shinydashboard","ggplot2","dplyr","tidyr","DT","ggrepel"))

library(shiny)
library(shinydashboard)
library(ggplot2)
library(dplyr)
library(tidyr)
library(DT)

# ============================================================
# LOAD AND PREPARE REAL GDSC DATA
# ============================================================

cat("Loading real GDSC2 data...\n")

gdsc_raw <- read.csv("GDSC2_fitted_dose_response_27Oct23.csv",
                     stringsAsFactors = FALSE)

# Filter to well-annotated cancer types
gdsc <- gdsc_raw %>%
  filter(
    !is.na(LN_IC50),
    !is.na(AUC),
    TCGA_DESC != "UNCLASSIFIED",
    WEBRELEASE == "Y"
  ) %>%
  mutate(
    IC50       = exp(LN_IC50),
    cancer_type = TCGA_DESC,
    drug        = DRUG_NAME,
    cell_line   = CELL_LINE_NAME,
    pathway     = PATHWAY_NAME,
    target      = PUTATIVE_TARGET
  )

# Top cancer types and drugs for UI
cancer_counts <- as.data.frame(table(gdsc$cancer_type))
cancer_counts <- cancer_counts[order(-cancer_counts$Freq), ]
top_cancers   <- as.character(cancer_counts$Var1[1:15])

drug_counts <- as.data.frame(table(gdsc$drug))
drug_counts <- drug_counts[order(-drug_counts$Freq), ]
top_drugs   <- as.character(drug_counts$Var1[1:50])

top_pathways <- sort(unique(gdsc$pathway))

cat("Data loaded:", nrow(gdsc), "measurements\n")
cat("Cancer types:", length(unique(gdsc$cancer_type)), "\n")
cat("Drugs:", length(unique(gdsc$drug)), "\n")

# ============================================================
# UI
# ============================================================

ui <- dashboardPage(
  skin = "blue",

  dashboardHeader(title = "Drug Sensitivity Explorer — Real GDSC2 Data"),

  dashboardSidebar(
    sidebarMenu(
      menuItem("Drug Explorer",    tabName = "drug",    icon = icon("pills")),
      menuItem("Cancer Analysis",  tabName = "cancer",  icon = icon("microscope")),
      menuItem("Pathway Analysis", tabName = "pathway", icon = icon("project-diagram")),
      menuItem("Data Table",       tabName = "table",   icon = icon("table"))
    ),
    hr(),
    p(style = "padding-left:15px; color:#aaa; font-size:11px;",
      "Source: GDSC2 Release 8.5", br(),
      "Sanger Institute", br(),
      "286 drugs | 969 cell lines")
  ),

  dashboardBody(
    tags$head(tags$style(HTML("
      .box { border-radius: 6px; }
      .value-box { border-radius: 6px; }
      .small-box { border-radius: 6px; }
    "))),

    tabItems(

      # ── TAB 1: DRUG EXPLORER ──
      tabItem(tabName = "drug",
        fluidRow(
          box(width = 3, title = "Controls", status = "primary", solidHeader = TRUE,
            selectInput("drug_select", "Select Drug:",
                        choices  = sort(top_drugs),
                        selected = "Erlotinib"),
            hr(),
            p(style="font-size:12px; color:#666;",
              "Select a drug to explore its sensitivity profile across cancer types.
               Lower IC50 = more sensitive.")
          ),
          box(width = 9, title = "Drug Sensitivity by Cancer Type", status = "primary", solidHeader = TRUE,
            plotOutput("drug_boxplot", height = "420px")
          )
        ),
        fluidRow(
          box(width = 6, title = "IC50 Distribution", status = "info", solidHeader = TRUE,
            plotOutput("drug_hist", height = "280px")
          ),
          box(width = 6, title = "Drug Information", status = "info", solidHeader = TRUE,
            tableOutput("drug_info")
          )
        )
      ),

      # ── TAB 2: CANCER ANALYSIS ──
      tabItem(tabName = "cancer",
        fluidRow(
          box(width = 3, title = "Controls", status = "success", solidHeader = TRUE,
            selectInput("cancer_select", "Select Cancer Type:",
                        choices  = sort(top_cancers),
                        selected = "LUAD"),
            selectInput("metric_select", "Sensitivity Metric:",
                        choices  = c("LN_IC50", "AUC", "Z_SCORE"),
                        selected = "LN_IC50"),
            sliderInput("top_n", "Number of top drugs:",
                        min = 5, max = 30, value = 15, step = 5),
            hr(),
            p(style="font-size:12px; color:#666;",
              "Lower LN_IC50 = more sensitive.", br(),
              "Lower AUC = more sensitive.", br(),
              "Lower Z_SCORE = more sensitive.")
          ),
          box(width = 9, title = "Most Sensitive Drugs", status = "success", solidHeader = TRUE,
            plotOutput("cancer_drugs", height = "450px")
          )
        ),
        fluidRow(
          box(width = 12, title = "Sensitivity Heatmap — Top Drugs vs Cancer Types",
              status = "success", solidHeader = TRUE,
            plotOutput("heatmap", height = "380px")
          )
        )
      ),

      # ── TAB 3: PATHWAY ANALYSIS ──
      tabItem(tabName = "pathway",
        fluidRow(
          box(width = 3, title = "Controls", status = "warning", solidHeader = TRUE,
            selectInput("pathway_select", "Select Pathway:",
                        choices  = top_pathways,
                        selected = "EGFR signaling"),
            selectInput("cancer_pathway", "Cancer Type:",
                        choices  = c("All", sort(top_cancers)),
                        selected = "All"),
            hr(),
            p(style="font-size:12px; color:#666;",
              "Explore drug sensitivity by targeted pathway.")
          ),
          box(width = 9, title = "Drug Sensitivity by Pathway", status = "warning", solidHeader = TRUE,
            plotOutput("pathway_plot", height = "420px")
          )
        ),
        fluidRow(
          box(width = 6, title = "Pathway Overview", status = "warning", solidHeader = TRUE,
            plotOutput("pathway_overview", height = "300px")
          ),
          box(width = 6, title = "Pathway Statistics", status = "warning", solidHeader = TRUE,
            tableOutput("pathway_stats")
          )
        )
      ),

      # ── TAB 4: DATA TABLE ──
      tabItem(tabName = "table",
        fluidRow(
          box(width = 12, title = "GDSC2 Data Explorer", status = "primary", solidHeader = TRUE,
            fluidRow(
              column(3, selectInput("table_cancer", "Cancer Type:",
                                   choices = c("All", sort(top_cancers)),
                                   selected = "LUAD")),
              column(3, selectInput("table_drug", "Drug:",
                                   choices = c("All", sort(top_drugs)),
                                   selected = "All")),
              column(3, selectInput("table_pathway", "Pathway:",
                                   choices = c("All", top_pathways),
                                   selected = "All")),
              column(3, br(), downloadButton("download_data", "Download CSV",
                                             class = "btn-primary"))
            ),
            hr(),
            DTOutput("data_table")
          )
        )
      )
    )
  )
)

# ============================================================
# SERVER
# ============================================================

server <- function(input, output, session) {

  # ── Reactive: filter by drug ──
  drug_data <- reactive({
    gdsc %>% filter(drug == input$drug_select)
  })

  # ── Reactive: filter by cancer ──
  cancer_data <- reactive({
    gdsc %>% filter(cancer_type == input$cancer_select)
  })

  # ── Reactive: filter by pathway ──
  pathway_data <- reactive({
    df <- gdsc %>% filter(pathway == input$pathway_select)
    if (input$cancer_pathway != "All")
      df <- df %>% filter(cancer_type == input$cancer_pathway)
    df
  })

  # ── Tab 1: Drug boxplot ──
  output$drug_boxplot <- renderPlot({
    df <- drug_data()
    if (nrow(df) == 0) return(NULL)

    df_top <- df %>%
      group_by(cancer_type) %>%
      filter(n() >= 3) %>%
      ungroup()

    med_order <- df_top %>%
      group_by(cancer_type) %>%
      summarise(med = median(LN_IC50, na.rm=TRUE)) %>%
      arrange(med)

    df_top$cancer_type <- factor(df_top$cancer_type,
                                  levels = med_order$cancer_type)

    ggplot(df_top, aes(x=cancer_type, y=LN_IC50, fill=cancer_type)) +
      geom_boxplot(alpha=0.8, outlier.size=0.8) +
      geom_hline(yintercept=0, linetype="dashed", color="grey50") +
      coord_flip() +
      scale_fill_viridis_d(option="turbo", alpha=0.8) +
      labs(
        title   = paste(input$drug_select, "— Sensitivity across Cancer Types"),
        subtitle = paste("Real GDSC2 data |", nrow(df_top), "cell lines | Lower LN_IC50 = more sensitive"),
        x = "Cancer Type", y = "LN_IC50",
        caption = "Source: GDSC2 Release 8.5, Sanger Institute"
      ) +
      theme_classic(base_size=12) +
      theme(legend.position="none",
            plot.title=element_text(face="bold"))
  })

  # ── Tab 1: Drug histogram ──
  output$drug_hist <- renderPlot({
    df <- drug_data()
    if (nrow(df) == 0) return(NULL)
    ggplot(df, aes(x=LN_IC50)) +
      geom_histogram(bins=40, fill="#2A9D8F", alpha=0.8, color="white") +
      geom_vline(xintercept=median(df$LN_IC50, na.rm=TRUE),
                 color="#E76F51", linewidth=1.2, linetype="dashed") +
      labs(title=paste(input$drug_select, "— IC50 Distribution"),
           x="LN_IC50", y="Count",
           subtitle=paste("Median LN_IC50 =",
                          round(median(df$LN_IC50, na.rm=TRUE), 2))) +
      theme_classic(base_size=12) +
      theme(plot.title=element_text(face="bold"))
  })

  # ── Tab 1: Drug info table ──
  output$drug_info <- renderTable({
    df <- drug_data()
    if (nrow(df) == 0) return(NULL)
    data.frame(
      Metric = c("Drug", "Target", "Pathway", "Cell lines tested",
                 "Cancer types", "Median LN_IC50", "Median AUC"),
      Value  = c(
        unique(df$drug)[1],
        unique(df$target)[1],
        unique(df$pathway)[1],
        as.character(nrow(df)),
        as.character(length(unique(df$cancer_type))),
        round(median(df$LN_IC50, na.rm=TRUE), 3),
        round(median(df$AUC, na.rm=TRUE), 3)
      )
    )
  })

  # ── Tab 2: Cancer drugs barplot ──
  output$cancer_drugs <- renderPlot({
    df <- cancer_data()
    if (nrow(df) == 0) return(NULL)

    metric <- input$metric_select
    top_df <- df %>%
      group_by(drug, pathway) %>%
      summarise(mean_val = mean(.data[[metric]], na.rm=TRUE),
                n_lines  = n(), .groups="drop") %>%
      arrange(mean_val) %>%
      as.data.frame()
    top_df <- top_df[1:min(input$top_n, nrow(top_df)), ]

    ggplot(top_df, aes(x=reorder(drug, mean_val), y=mean_val, fill=pathway)) +
      geom_bar(stat="identity", alpha=0.85) +
      coord_flip() +
      scale_fill_brewer(palette="Set3") +
      labs(
        title   = paste("Most Sensitive Drugs in", input$cancer_select),
        subtitle = paste("Real GDSC2 data | Metric:", metric,
                         "| Lower = more sensitive"),
        x = "Drug", y = metric, fill = "Pathway",
        caption = "Source: GDSC2 Release 8.5, Sanger Institute"
      ) +
      theme_classic(base_size=12) +
      theme(plot.title=element_text(face="bold"),
            legend.position="bottom",
            legend.text=element_text(size=8))
  })

  # ── Tab 2: Heatmap ──
  output$heatmap <- renderPlot({
    # Top 10 drugs and top 10 cancer types using base R
    drug_freq     <- sort(table(gdsc$drug), decreasing=TRUE)
    top10_drugs   <- names(drug_freq)[1:10]
    cancer_freq   <- sort(table(gdsc$cancer_type[gdsc$cancer_type %in% top_cancers]),
                          decreasing=TRUE)
    top10_cancers <- names(cancer_freq)[1:10]

    heat_df <- gdsc %>%
      filter(drug %in% top10_drugs, cancer_type %in% top10_cancers) %>%
      group_by(drug, cancer_type) %>%
      summarise(mean_ic50 = mean(LN_IC50, na.rm=TRUE), .groups="drop")

    ggplot(heat_df, aes(x=cancer_type, y=drug, fill=mean_ic50)) +
      geom_tile(color="white", linewidth=0.4) +
      scale_fill_gradient2(low="#2A9D8F", mid="white", high="#E76F51",
                            midpoint=median(heat_df$mean_ic50, na.rm=TRUE),
                            name="Mean\nLN_IC50") +
      labs(
        title   = "Drug Sensitivity Heatmap — Top Drugs vs Cancer Types",
        subtitle = "Real GDSC2 data | Lower LN_IC50 = more sensitive (green)",
        x = "Cancer Type", y = "Drug",
        caption = "Source: GDSC2 Release 8.5, Sanger Institute"
      ) +
      theme_classic(base_size=11) +
      theme(plot.title=element_text(face="bold"),
            axis.text.x=element_text(angle=45, hjust=1))
  })

  # ── Tab 3: Pathway plot ──
  output$pathway_plot <- renderPlot({
    df <- pathway_data()
    if (nrow(df) == 0) return(NULL)

    # Limit to top 20 most tested drugs to avoid overcrowding
    top_drugs_pathway <- df %>%
      group_by(drug) %>%
      summarise(n = n(), mean_ic50 = mean(LN_IC50, na.rm=TRUE), .groups="drop") %>%
      arrange(mean_ic50) %>%
      as.data.frame()
    top_drugs_pathway <- top_drugs_pathway[1:min(20, nrow(top_drugs_pathway)), ]

    df_filtered <- df %>%
      filter(drug %in% top_drugs_pathway$drug)

    n_drugs <- length(unique(df_filtered$drug))
    label_size <- ifelse(n_drugs > 15, 8, ifelse(n_drugs > 10, 9, 10))

    ggplot(df_filtered, aes(x=reorder(drug, LN_IC50, FUN=median),
                             y=LN_IC50, fill=cancer_type)) +
      geom_boxplot(alpha=0.75, outlier.size=0.5) +
      coord_flip() +
      scale_fill_brewer(palette="Set2") +
      labs(
        title    = paste("Pathway:", input$pathway_select),
        subtitle = paste("Top 20 most sensitive drugs | Real GDSC2 data"),
        x="Drug", y="LN_IC50", fill="Cancer Type",
        caption="Source: GDSC2 Release 8.5, Sanger Institute"
      ) +
      theme_classic(base_size=12) +
      theme(plot.title=element_text(face="bold"),
            axis.text.y=element_text(size=label_size),
            legend.position="bottom",
            legend.text=element_text(size=8))
  })

  # ── Tab 3: Pathway overview ──
  output$pathway_overview <- renderPlot({
    pathway_counts <- gdsc %>%
      group_by(pathway) %>%
      summarise(n_drugs=n_distinct(drug),
                mean_ic50=mean(LN_IC50, na.rm=TRUE),
                .groups="drop") %>%
      arrange(desc(n_drugs)) %>%
      as.data.frame()
    pathway_counts <- pathway_counts[1:min(12, nrow(pathway_counts)), ]

    ggplot(pathway_counts, aes(x=reorder(pathway, n_drugs),
                                y=n_drugs, fill=mean_ic50)) +
      geom_bar(stat="identity", alpha=0.85) +
      coord_flip() +
      scale_fill_gradient(low="#2A9D8F", high="#E76F51", name="Mean\nLN_IC50") +
      labs(title="Drugs per Pathway",
           x="Pathway", y="Number of Drugs") +
      theme_classic(base_size=11) +
      theme(plot.title=element_text(face="bold"))
  })

  # ── Tab 3: Pathway stats ──
  output$pathway_stats <- renderTable({
    df <- pathway_data()
    if (nrow(df) == 0) return(NULL)
    data.frame(
      Metric = c("Pathway", "Drugs tested", "Cell lines",
                 "Cancer types", "Mean LN_IC50", "Mean AUC"),
      Value  = c(
        input$pathway_select,
        as.character(length(unique(df$drug))),
        as.character(nrow(df)),
        as.character(length(unique(df$cancer_type))),
        round(mean(df$LN_IC50, na.rm=TRUE), 3),
        round(mean(df$AUC, na.rm=TRUE), 3)
      )
    )
  })

  # ── Tab 4: Data table ──
  filtered_table <- reactive({
    df <- gdsc %>%
      select(cancer_type, cell_line, drug, pathway, target,
             LN_IC50, AUC, Z_SCORE) %>%
      mutate(across(where(is.numeric), ~round(., 3)))
    if (input$table_cancer != "All") df <- df %>% filter(cancer_type==input$table_cancer)
    if (input$table_drug   != "All") df <- df %>% filter(drug==input$table_drug)
    if (input$table_pathway != "All") df <- df %>% filter(pathway==input$table_pathway)
    df
  })

  output$data_table <- renderDT({
    datatable(filtered_table(),
              options=list(pageLength=15, scrollX=TRUE),
              filter="top", rownames=FALSE)
  })

  output$download_data <- downloadHandler(
    filename = function() paste0("GDSC2_", Sys.Date(), ".csv"),
    content  = function(file) write.csv(filtered_table(), file, row.names=FALSE)
  )
}

shinyApp(ui, server)
