library(data.table)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(scales)
library(viridis)
library(gridExtra)
library(patchwork)
library(RColorBrewer)
library(ggpmisc)
library(ggrastr)
library(BAS)
library(R.utils)
theme_set(theme_cowplot(font_size = 9))


#########################
### Set the following to the repo directory
setwd('~/research/specificity_length_luck/')


#########################
### Figures 1A and 1B ###
#########################
set.seed(42)

N=500
B=1

df=data.frame(g=rbinom(N,2,0.15))
df$y=B*df$g+rnorm(N,0,1)
df$D=0
df[df$y>1,]$D=1
df$g=factor(df$g,levels = c(0,1,2))

mod1=lm(y~g,df)

p1<-ggplot(df, aes(x=g, y=y)) +
  #geom_boxplot(fill='#BADBFF', color="black",outlier.size = 0.3)+
  #geom_point(color="black",size = 1)+
  geom_jitter(width = 0.1,size = 0.5, color=alpha("#FF5959", 0.3))+
  geom_abline( slope= 1,intercept = -1, linetype="solid", linewidth=1)+
  xlab("Genotype") +
  ylab("Phenotype") +
  scale_x_discrete(labels = c("AA", "AC", "CC")) +
  theme_cowplot(9)+
  theme(axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

outfile="paper_figures/fig1a.pdf"
save_plot(outfile, p1,base_width = 2.3, base_height = 2)


#######

N=500
B=-1

df=data.frame(g=rbinom(N,1,0.05))
df$y=B*df$g+rnorm(N,0,1)
df$g=factor(df$g,levels = c(0,1))

mod1=lm(y~g,df)

p1<-ggplot(df, aes(x=g, y=y)) +
  #geom_boxplot(fill='#BADBFF', color="black",outlier.size = 0.3)+
  #geom_point(color="black",size = 1)+
  geom_jitter(width = 0.0666,size = 0.5,color=alpha("#1E56A0", 0.3)) +
  geom_abline( slope= -1,intercept = 1, linetype="solid", linewidth=1)+
  #xlab("") +
  ylab("Phenotype") +
  theme_cowplot(9)+
  scale_x_discrete(labels = c("individuals with\n no LoFs in gene", "individuals with\n any LoF in gene")) +
  theme(axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.y = element_blank())

outfile="paper_figures/fig1b.pdf"
save_plot(outfile, p1,base_width = 2.3, base_height = 2)











#########################
### Figures 1C,D, S1-3 ###
#########################

d_all_gwas_hits=fread("data/all_gwas_hits.txt")
d_all_top_hits=fread("data/top_gwas_hits.txt")
coloc_gene_counts=fread("data/coloc_status.txt")

#####
#select traits to plot

d_trait_names=fread("data/trait_lookup.tsv",header = F)
colnames(d_trait_names)=c("pheno","trait")

d_indep_traits=fread("data/indep_traits.txt")
colnames(d_indep_traits)="pheno"
d_indep_traits=left_join(d_indep_traits,d_trait_names)

d_cont_traits=data.frame(pheno=unique(d_all_gwas_hits$pheno))
d_cont_traits=left_join(d_cont_traits,d_trait_names)

coloc_gene_counts=left_join(coloc_gene_counts,d_cont_traits)

coloc_gene_counts[coloc_gene_counts$pheno %in% d_indep_traits$pheno,]
coloc_gene_counts[coloc_gene_counts$burden_genes>=10,] %>% select(pheno,trait,burden_genes)
coloc_gene_counts[coloc_gene_counts$burden_genes>=10 & (coloc_gene_counts$pheno %in% d_indep_traits$pheno),]
coloc_gene_counts[coloc_gene_counts$burden_genes>=10,] %>% arrange(trait) %>% select(trait,burden_genes)

#######

selected_traits=c("50_irnt","3062_irnt","23106_irnt",
                  "30780_irnt","30760_irnt","23105_irnt",
                  "30880_irnt","30830_irnt","30770_irnt",
                  "30050_irnt","30100_irnt","30010_irnt")

coloc_selected=coloc_gene_counts[coloc_gene_counts$pheno %in% selected_traits,]
coloc_selected$name=c("FVC","LDL","WBI","SHBG","IGF-1","Urate","BMR","HDL","RBC count","MPV","MCH","Height")

coloc_selected$name=factor(coloc_selected$name,levels = rev(coloc_selected$name))
coloc_selected$type <- factor(coloc_selected$name, levels = c("Burden", "Coloc"))

p1<-ggplot(coloc_selected, aes(x = name, y = burden_genes, fill = "Burden only")) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_bar(aes(y = coloc_count, fill = "Burden and GWAS"), stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c("Burden and GWAS" = "#1f77b4", "Burden only" = "#aec7e8"), breaks = c("", "")) +
  # scale_fill_manual(values = c("Burden and GWAS" = "#1f77b4", "Burden only" = "#aec7e8"),
  #                   breaks = c("Burden only", "Burden and GWAS")) +
  labs(x = "", y = "Number of genes") +
  theme_cowplot(9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        legend.position = c(0.3, 0.9),
        legend.title = element_blank())

outfile=paste0("paper_figures/fig1c.pdf")
save_plot(outfile, p1,base_width = 3, base_height = 3)


#######
#scatter plots
d_z_scores=fread("data/zscores_by_trait_MAF1_hgncIDs.txt")

spearmans <- list()

for (trait_temp in unique(d_all_gwas_hits$pheno)){
  
  d_burden_temp=d_z_scores %>% select(trait_temp)
  colnames(d_burden_temp)="z"
  gene_num=nrow(d_burden_temp[is.finite(d_burden_temp$z) & !is.na(d_burden_temp$z),])
  
  df=d_all_gwas_hits[d_all_gwas_hits$pheno==trait_temp,]
  df$burden_pval=2 * pnorm(abs(df$max_burden_z), lower.tail = FALSE)
  
  #df[df$pval==0,]$pval=min(1e-300,min(df[df$pval>0]$pval,na.rm = T))
  #df[df$burden_pval==0,]$burden_pval=min(1e-300,min(df[df$burden_pval>0]$burden_pval,na.rm = T))
  
  df[df$pval==0,]$pval=min(df[df$pval>0]$pval,na.rm = T)
  df[df$burden_pval==0,]$burden_pval=min(df[df$burden_pval>0]$burden_pval,na.rm = T)
  
  df$log_gwas_p=-log10(df$pval)
  df$log_burden_p = -log10(df$burden_pval)
  
  burden_cutoff=-log10(0.05/gene_num)
  gwas_cutoff=-log10(5e-8)
  
  if(trait_temp == "50_irnt"){
    p1 <- ggplot(df, aes(x=log_gwas_p, y=log_burden_p)) +
      geom_point(size=1.1,color = "black",alpha=0.65) +
      geom_hline(yintercept = burden_cutoff, color = "red", linetype = "dashed") +
      geom_vline(xintercept = gwas_cutoff, linetype = "dashed",color = "red") +
      scale_x_continuous(breaks = seq(0, 302, by = 50)) +
      xlab(expression(paste(-log[10], " p (GWAS)"))) +
      ylab(expression(paste(-log[10], " min p (burden) in GWAS locus"))) +
      theme_cowplot(9)
    outfile=paste0("paper_figures/fig1d.pdf")
    save_plot(outfile, p1,base_width = 3, base_height = 3)
  }
  spearmans[[trait_temp]] <- cor(df$log_gwas_p, df$log_burden_p, method='spearman')
}



#### supp plot of spearman's across traits
spearman_hist <- ggplot(data.frame(unlist(spearmans)), aes(x=unlist.spearmans.)) +
  geom_histogram(bins=60) +
  xlab(expression(paste("Spearman's correlation between ", -log[10], ' p (GWAS) and ',-log[10], ' min p (burden) in GWAS locus'))) + 
  ylab("Number of traits") +
  theme_cowplot(9)

outfile=paste0("paper_figures/figS3.pdf")
save_plot(outfile, spearman_hist, base_width = 8, base_height = 3)


#### supp plot with many traits in grid
to_plot <- list()
for(idx in 1:10){
  trait_temp <- c("30050_irnt",
                  "30100_irnt",
                  "30010_irnt",
                  "30760_irnt",
                  "23105_irnt",
                  "30880_irnt",
                  "30770_irnt",
                  "30830_irnt",
                  "23106_irnt",
                  "30780_irnt")[idx]
  trait_name <- c("MCH",
                  "MPV",
                  "RBC count",
                  "HDL",
                  "BMR",
                  "Urate",
                  "IGF-1",
                  "SHBG",
                  "WBI",
                  "LDL")[idx]
  d_burden_temp=d_z_scores %>% select(trait_temp)
  colnames(d_burden_temp)="z"
  gene_num=nrow(d_burden_temp[is.finite(d_burden_temp$z) & !is.na(d_burden_temp$z),])
  
  df=d_all_gwas_hits[d_all_gwas_hits$pheno==trait_temp,]
  df$burden_pval=2 * pnorm(abs(df$max_burden_z), lower.tail = FALSE)
  
  #df[df$pval==0,]$pval=min(1e-300,min(df[df$pval>0]$pval,na.rm = T))
  #df[df$burden_pval==0,]$burden_pval=min(1e-300,min(df[df$burden_pval>0]$burden_pval,na.rm = T))
  
  df[df$pval==0,]$pval=min(df[df$pval>0]$pval,na.rm = T)
  df[df$burden_pval==0,]$burden_pval=min(df[df$burden_pval>0]$burden_pval,na.rm = T)
  
  df$log_gwas_p=-log10(df$pval)
  df$log_burden_p = -log10(df$burden_pval)
  
  burden_cutoff=-log10(0.05/gene_num)
  gwas_cutoff=-log10(5e-8)
  
  p1 <- ggplot(df, aes(x=log_gwas_p, y=log_burden_p)) +
    geom_point(size=1.1,color = "black",alpha=0.65) +
    geom_hline(yintercept = burden_cutoff, color = "red", linetype = "dashed") +
    geom_vline(xintercept = gwas_cutoff, linetype = "dashed",color = "red") +
    scale_x_continuous(breaks = seq(0, 302, by = 50)) +
    xlab(expression(paste(-log[10], " p (GWAS)"))) +
    ylab(expression(paste(-log[10], " min p (burden) in GWAS locus"))) +
    theme_cowplot(9) +
    ggtitle(trait_name)
  
  to_plot[[idx]] <- p1
}

to_plot_full <- plot_grid(to_plot[[1]], to_plot[[2]], to_plot[[3]], to_plot[[4]], to_plot[[5]],
                          to_plot[[6]], to_plot[[7]], to_plot[[8]], to_plot[[9]],
                          nrow=3, ncol=3, byrow=TRUE)

save_plot('paper_figures/figS2.pdf', to_plot_full, base_width = 3*3, base_height = 3*3)



coloc_gene_counts$name = paste(coloc_gene_counts$pheno,  coloc_gene_counts$trait, sep=': ')
coloc_gene_counts$name=factor(coloc_gene_counts$name,levels = coloc_gene_counts$name)
coloc_gene_counts$type <- factor(coloc_gene_counts$name, levels = c("Burden", "Coloc"))

p1_supp<-ggplot(coloc_gene_counts, aes(x = name, y = burden_genes, fill = "Burden only")) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_bar(aes(y = coloc_count, fill = "Burden and GWAS"), stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c("Burden and GWAS" = "#1f77b4", "Burden only" = "#aec7e8"), breaks = c("", "")) +
  # scale_fill_manual(values = c("Burden and GWAS" = "#1f77b4", "Burden only" = "#aec7e8"),
  #                   breaks = c("Burden only", "Burden and GWAS")) +
  labs(x = "", y = "Number of genes") +
  theme_cowplot(7) +
  coord_flip() + 
  theme(# axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.3, 0.9),
    legend.title = element_blank())

outfile=paste0("paper_figures/figS1.pdf")
save_plot(outfile, p1_supp, base_width = 2.25*4, base_height = 2*5)










#########################
### Figure 1E ###
#########################

d_blocks=fread("data/fourier_ls-all.named.bed")
colnames(d_blocks)=c("chr","start","end","block")
d_blocks$chr=sub("^chr", "", d_blocks$chr)
d_blocks$chr=as.numeric(d_blocks$chr)

d_gene_blocks=fread("data/genes.protein_coding.blocks.v39.txt",header = F)
colnames(d_gene_blocks)=c("gene","tss_block","tes_block")
bad_blocks1=d_gene_blocks[d_gene_blocks$tss_block!=d_gene_blocks$tes_block,]$tss_block
bad_blocks2=d_gene_blocks[d_gene_blocks$tss_block!=d_gene_blocks$tes_block,]$tes_block
bad_blocks=unique(c(bad_blocks1,bad_blocks2))
d_gene_blocks$gene_block=d_gene_blocks$tss_block
d_gene_blocks$tss_tes_block=paste0(d_gene_blocks$tss_block,"_",d_gene_blocks$tes_block)
d_gene_blocks[d_gene_blocks$tss_block!=d_gene_blocks$tes_block,]$gene_block=d_gene_blocks[d_gene_blocks$tss_block!=d_gene_blocks$tes_block,]$tss_tes_block
d_gene_blocks=d_gene_blocks %>% select(gene,gene_block)

d_hits=fread("data/gwas_hits.blocks.txt")

d_burden_z=fread("data/zscores_by_trait_MAF1_hgncIDs.txt")
colnames(d_burden_z)[1]="gene"
#d_burden_z=left_join(d_gene_blocks,d_burden_z)


##########
trait="50_irnt"

d_gwas=d_hits[d_hits$pheno==trait,]
d_gwas[d_gwas$pval==0,]$pval=1e-300
#d_gwas=left_join(d_gwas,d_gene_blocks)

d_burden=d_burden_z %>% select(gene,trait) 
colnames(d_burden)[2]="trait"
d_burden= d_burden%>% arrange(-abs(trait))


##

d_genes=fread("data/genes.protein_coding.v39.gtf")

dx=d_genes %>% select(chr,hgnc_id)
colnames(dx)=c("chr_str","gene")

dy=d_genes %>% select(chr,gene,hgnc_id,start,end)
colnames(dy)=c("chr_str","gene_name","gene","start","end")


window_size = 100000

head(d_burden)

gwas_snp="9:35932874:G:A"
gene_temp="HGNC:7944"
gwas_chr="chr9"
gwas_pos=floor((d_genes[d_genes$hgnc_id==gene_temp,]$start + d_genes[d_genes$hgnc_id==gene_temp,]$end)/2)
block_temp="block_963"

d_gwas_signal=left_join(d_gwas,dx)
d_gwas_signal=d_gwas_signal %>% separate(SNP,into = c("chr","pos","a1","a2"),sep = ":",remove = F)
d_gwas_signal$pos=as.numeric(d_gwas_signal$pos)
d_gwas_signal=d_gwas_signal[d_gwas_signal$chr_str==gwas_chr,]
d_gwas_signal=d_gwas_signal[d_gwas_signal$pos<(gwas_pos + window_size) & d_gwas_signal$pos>(gwas_pos - window_size),]

d_burden_signal=left_join(d_burden,dy)
d_burden_signal=d_burden_signal[d_burden_signal$chr_str==gwas_chr,]
d_burden_signal1=d_burden_signal[d_burden_signal$start<(gwas_pos + window_size) & d_burden_signal$start>(gwas_pos - window_size),]
d_burden_signal2=d_burden_signal[d_burden_signal$end<(gwas_pos + window_size) & d_burden_signal$end>(gwas_pos - window_size),]
d_burden_signal=d_burden_signal[d_burden_signal$gene %in% intersect(d_burden_signal1$gene,d_burden_signal2$gene),]
d_burden_signal$pval <- 2 * pnorm(-abs(d_burden_signal$trait))
d_burden_signal[is.na(d_burden_signal$trait),]$trait=1
d_burden_signal$abs_z=abs(d_burden_signal$trait)

max_burden_signal= min(2 * pnorm(-abs(d_burden$trait)),na.rm = T)
max_burden_signal= max(abs(d_burden$trait),na.rm = T)

z_scores_cutoff <- qnorm(1 - 0.05/nrow(d_burden)/2) 


sequence_legend=c(0, 5, 10, 15, 20, 25)
grad_labels = as.character(sequence_legend)
grad_labels[6]=">25"

#


loc_genes=d_burden_signal$gene
d_gene_loc=dy[dy$gene %in% loc_genes] %>% select(start,end,gene_name)
d_gene_loc$start=d_gene_loc$start/1000000
d_gene_loc$end=d_gene_loc$end/1000000

##########


# Sample data
x <- (d_gwas_signal$pos)/1000000  # Genomic coordinates
p <- -log10(d_gwas_signal$pval)    # Signal for experiment 1

d_temp=left_join(d_gene_loc,d_burden_signal %>% select(abs_z,gene_name))
q <- -log10(2*pnorm(-d_temp$abs_z))    # Signal for experiment 2
q[is.na(q)] <- 0
q[q>25] <- 25
# q[q<4.7] <- 0



# Function to calculate stacking order based on padded start and end positions
calculate_stacking_order <- function(data, window_size) {
  data$start <- data$start - window_size  # Pad start positions
  data$end <- data$end + window_size  # Pad end positions
  
  stack_order <- rep(1, nrow(data))  # Assign all genes to stack 1
  
  for (k in 2:nrow(data)) {
    for (j in 1:(k - 1)) {
      # Check for overlap within the same stack
      if (stack_order[k] == stack_order[j] && 
          data$start[k] < data$end[j] && data$end[k] > data$start[j]) {
        stack_order[k] <- stack_order[k] + 1  # Increment stack number
        j <- 0  # Reset j to start from the beginning
      }
    }
  }
  
  return(stack_order)  # Return stack order
}


# Sample gene locations
gene_data <- d_gene_loc
gene_data$stack_order <- calculate_stacking_order(gene_data,0.02)
gene_data$bold_indicator=0
gene_data[gene_data$gene_name=="NPR2",]$bold_indicator=1

# Create a dataframe for ggplot
df <- data.frame(
  x = x,
  p = p
)

# Define color gradient function for genes
color_gradient <- function(value) {
  color_intensity <- value/10  # Reverse intensity: small q values correspond to white
  color_intensity <- color_intensity * 255  # Convert to RGB scale
  rgb(color_intensity, 0, 0, maxColorValue = 255)
}

max_x <- max(max(df$x), max(gene_data$end))
min_x <- min(min(df$x), min(gene_data$start))

#########


# Define gene_plot with desired changes
gene_plot <- ggplot(gene_data) +
  geom_rect(aes(xmin = start, xmax = end, ymin = -stack_order / 10 - 0.015, ymax = -stack_order / 10 + 0.015, fill = q), color = "black") +  # Gene bodies as rectangles with black borders
  geom_text(aes(x = (start + end) / 2, y = -stack_order / 10 - 0.02, label = gene_name, group = stack_order, fontface = ifelse(bold_indicator == 1, "bold", "plain")), size = 2, color = "black", vjust = 1.2, hjust = 0.5) +  # Gene labels with adjusted vertical positioning and black color, stacked systematically
  scale_fill_gradient2(low = "white", high = "#FF5733", breaks = sequence_legend,labels=grad_labels, limits = range(sequence_legend))+       ## color of the corresponding aes
  coord_cartesian(xlim = c(min_x-0.02, max_x+0.02)) +
  theme_void()+
  theme(plot.background = element_blank(),  # Set plot background to white
        panel.grid = element_blank(),  # Remove grid lines
        axis.title.x = element_blank(),  # Remove x-axis title
        axis.title.y = element_blank(),  # Remove y-axis title
        axis.text.y = element_blank(),  # Remove y-axis text
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        axis.line = element_blank(),  # Remove axis lines
        legend.position = "left",  # Move legend to the right
        # plot.margin = unit(rep(1, 4),'lines'),
        legend.title = element_text(size = 9, angle = 90, hjust=0.6),   # Adjust legend title size
        legend.title.position = 'left',
        legend.text = element_text(size = 6)) +  # Adjust legend title size
  theme(legend.key.size = unit(0.5, "lines"), legend.position="left") +  # Adjust legend key size to make it smaller
  # labs(fill = "") + 
  labs(fill = expression(paste("  ", -log[10], " p (burden)"))) +  # Specify legend title with LaTeX formatting
  expand_limits(y = c(-0.56, -0.1))  # Adjust the plot margins

# Define main_plot
main_plot <- ggplot(df) +
  geom_point(aes(x = x, y = p), color = "blue",alpha=0.7) +  # Plot Signal 1 as points
  ylab(expression(paste(-log[10], " p ", "(GWAS)      "))) +
  xlab("Position on chr9 (Mb)") +
  coord_cartesian(xlim = c(min_x-0.02, max_x+0.02), ylim = c(0, 210)) +  # Adjust axis limits
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "black") +
  theme_cowplot(9) +
  theme(plot.background = element_blank(),  # Remove plot background
        panel.grid = element_blank(),  # Remove grid lines
        legend.position = "none") 

combined_plot <- (
  main_plot/gene_plot
) & theme(legend.position = c(-0.067, 0.4))

combined_plot
combined_plot <- combined_plot + plot_layout(heights = c(1, 2))

outfile="paper_figures/fig1e.pdf"
save_plot(outfile, combined_plot,base_width = 3, base_height = 2)












#########################
### Figure 1F ###
#########################

d_blocks=fread("data/fourier_ls-all.named.bed")
colnames(d_blocks)=c("chr","start","end","block")
d_blocks$chr=sub("^chr", "", d_blocks$chr)
d_blocks$chr=as.numeric(d_blocks$chr)

d_gene_blocks=fread("data/genes.protein_coding.blocks.v39.txt",header = F)
colnames(d_gene_blocks)=c("gene","tss_block","tes_block")
bad_blocks1=d_gene_blocks[d_gene_blocks$tss_block!=d_gene_blocks$tes_block,]$tss_block
bad_blocks2=d_gene_blocks[d_gene_blocks$tss_block!=d_gene_blocks$tes_block,]$tes_block
bad_blocks=unique(c(bad_blocks1,bad_blocks2))
d_gene_blocks$gene_block=d_gene_blocks$tss_block
d_gene_blocks$tss_tes_block=paste0(d_gene_blocks$tss_block,"_",d_gene_blocks$tes_block)
d_gene_blocks[d_gene_blocks$tss_block!=d_gene_blocks$tes_block,]$gene_block=d_gene_blocks[d_gene_blocks$tss_block!=d_gene_blocks$tes_block,]$tss_tes_block
d_gene_blocks=d_gene_blocks %>% select(gene,gene_block)

d_hits=fread("data/gwas_hits.blocks.txt")

d_burden_z=fread("data/zscores_by_trait_MAF1_hgncIDs.txt")
colnames(d_burden_z)[1]="gene"
#d_burden_z=left_join(d_gene_blocks,d_burden_z)


##########
trait="50_irnt"

d_gwas=d_hits[d_hits$pheno==trait,]
d_gwas[d_gwas$pval==0,]$pval=1e-300
#d_gwas=left_join(d_gwas,d_gene_blocks)

d_burden=d_burden_z %>% select(gene,trait) 
colnames(d_burden)[2]="trait"
d_burden= d_burden%>% arrange(-abs(trait))


##

d_genes=fread("data/genes.protein_coding.v39.gtf")

dx=d_genes %>% select(chr,hgnc_id)
colnames(dx)=c("chr_str","gene")

dy=d_genes %>% select(chr,gene,hgnc_id,start,end)
colnames(dy)=c("chr_str","gene_name","gene","start","end")


window_size = 1000000

head(d_gwas)

gwas_snp="4:145566477:T:C"
gene_temp="HGNC:14866"
gwas_chr="chr4"
gwas_pos=145566477
block_temp="block_493"

d_gwas_signal=left_join(d_gwas,dx)
d_gwas_signal=d_gwas_signal %>% separate(SNP,into = c("chr","pos","a1","a2"),sep = ":",remove = F)
d_gwas_signal$pos=as.numeric(d_gwas_signal$pos)
d_gwas_signal=d_gwas_signal[d_gwas_signal$chr_str==gwas_chr,]
d_gwas_signal=d_gwas_signal[d_gwas_signal$pos<(gwas_pos + window_size) & d_gwas_signal$pos>(gwas_pos - window_size),]
# d_gwas_signal[d_gwas_signal$pval==0,]$pval=1e-300

d_burden_signal=left_join(d_burden,dy)
d_burden_signal=d_burden_signal[d_burden_signal$chr_str==gwas_chr,]
d_burden_signal1=d_burden_signal[d_burden_signal$start<(gwas_pos + window_size) & d_burden_signal$start>(gwas_pos - window_size),]
d_burden_signal2=d_burden_signal[d_burden_signal$end<(gwas_pos + window_size) & d_burden_signal$end>(gwas_pos - window_size),]
d_burden_signal=d_burden_signal[d_burden_signal$gene %in% intersect(d_burden_signal1$gene,d_burden_signal2$gene),]
d_burden_signal$pval <- 2 * pnorm(-abs(d_burden_signal$trait))
d_burden_signal[is.na(d_burden_signal$trait),]$trait=1
d_burden_signal$abs_z=abs(d_burden_signal$trait)

max_burden_signal= min(2 * pnorm(-abs(d_burden$trait)),na.rm = T)
max_burden_signal= max(abs(d_burden$trait),na.rm = T)

z_scores_cutoff <- qnorm(1 - 0.05/nrow(d_burden)/2) 




sequence_legend=c(0, 5, 10, 15, 20, 25)
grad_labels = as.character(sequence_legend)
grad_labels[6]=">25"

#


loc_genes=d_burden_signal$gene
d_gene_loc=dy[dy$gene %in% loc_genes] %>% select(start,end,gene_name)
d_gene_loc$start=d_gene_loc$start/1000000
d_gene_loc$end=d_gene_loc$end/1000000

##########


# Sample data
x <- (d_gwas_signal$pos)/1000000  # Genomic coordinates
p <- -log10(d_gwas_signal$pval)    # Signal for experiment 1

d_temp=left_join(d_gene_loc,d_burden_signal %>% select(abs_z,gene_name))
q <- -log10(2*pnorm(-d_temp$abs_z))    # Signal for experiment 2
q[is.na(q)] <- 0
q[q>25] <- 25


# Function to calculate stacking order based on padded start and end positions
calculate_stacking_order <- function(data, window_size) {
  data$start <- data$start - window_size  # Pad start positions
  data$end <- data$end + window_size  # Pad end positions
  
  stack_order <- rep(1, nrow(data))  # Assign all genes to stack 1
  
  for (k in 2:nrow(data)) {
    for (j in 1:(k - 1)) {
      # Check for overlap within the same stack
      if (stack_order[k] == stack_order[j] && 
          data$start[k] < data$end[j] && data$end[k] > data$start[j]) {
        stack_order[k] <- stack_order[k] + 1  # Increment stack number
        j <- 0  # Reset j to start from the beginning
      }
    }
  }
  
  return(stack_order)  # Return stack order
}


# Sample gene locations
gene_data <- d_gene_loc

gene_data$stack_order=1
if (nrow(gene_data)>1){gene_data$stack_order <- calculate_stacking_order(gene_data,0.2)}

gene_data$bold_indicator=0
gene_data[gene_data$gene_name=="HHIP",]$bold_indicator=1

# Create a dataframe for ggplot
df <- data.frame(
  x = x,
  p = p
)

# Define color gradient function for genes
color_gradient <- function(value) {
  color_intensity <- 1 - value  # Reverse intensity: small q values correspond to white
  color_intensity <- color_intensity * 255  # Convert to RGB scale
  rgb(color_intensity, 0, 0, maxColorValue = 255)
}

max_x <- max(max(df$x), max(gene_data$end))
min_x <- min(min(df$x), min(gene_data$start))

#########


# Define gene_plot with desired changes
gene_plot <- ggplot(gene_data) +
  geom_rect(aes(xmin = start, xmax = end, ymin = -stack_order / 10 - 0.015, ymax = -stack_order / 10 + 0.015, fill = q), color = "black") +  # Gene bodies as rectangles with black borders
  geom_text(aes(x = (start + end) / 2, y = -stack_order / 10 - 0.02, label = gene_name, group = stack_order, fontface = ifelse(bold_indicator == 1, "bold", "plain")), size = 2, color = "black", vjust = 1.2, hjust = 0.5) +  # Gene labels with adjusted vertical positioning and black color, stacked systematically
  scale_fill_gradient2(low = "white", high = "#FF5733", breaks = sequence_legend,labels=grad_labels, limits = range(sequence_legend))+       ## color of the corresponding aes
  coord_cartesian(xlim = c(min_x-.02, max_x+.02)) +
  theme_void()+
  theme(plot.background = element_blank(),  # Set plot background to white
        panel.grid = element_blank(),  # Remove grid lines
        axis.title.x = element_blank(),  # Remove x-axis title
        axis.title.y = element_blank(),  # Remove y-axis title
        axis.text.y = element_blank(),  # Remove y-axis text
        axis.text.x = element_blank(),  # Remove x-axis text
        axis.ticks = element_blank(),  # Remove axis ticks
        axis.line = element_blank(),  # Remove axis lines
        legend.position = "left",  # Move legend to the right
        legend.title = element_text(size = 9, angle = 90, hjust=0.6),   # Adjust legend title size
        legend.title.position = 'left',
        legend.text = element_text(size = 6)) +  # Adjust legend title size
  theme(legend.key.size = unit(0.5, "lines")) +  # Adjust legend key size to make it smaller
  labs(fill = expression(paste("  ", -log[10], " p (burden)"))) +  # Specify legend title with LaTeX formatting
  expand_limits(y = c(-0.56, -0.1))  # Adjust the plot margins

# Define main_plot
main_plot <- ggplot(df) +
  geom_point(aes(x = x, y = p), color = "blue",alpha=0.7) +  # Plot Signal 1 as points
  ylab(expression(paste(-log[10], " p ", "(GWAS)      "))) +
  xlab("Position on chr4 (Mb)") +
  coord_cartesian(xlim = c(min_x-.02, max_x+.02), ylim = c(0, 210)) +  # Adjust axis limits
  geom_hline(yintercept = -log10(5e-8), linetype = "dashed", color = "black") +
  theme_cowplot(9) +
  theme(plot.background = element_blank(),  # Remove plot background
        panel.grid = element_blank(),  # Remove grid lines
        legend.position = "none") 

combined_plot <- (
  main_plot/gene_plot
) & theme(legend.position = c(-0.063, 0.4))

combined_plot <- combined_plot + plot_layout(heights = c(1, 2))

outfile="paper_figures/fig1f.pdf"
save_plot(outfile, combined_plot,base_width = 3, base_height = 2)






#########################
### Figure S4-7 ###
#########################


# scatter for height 
df <- read_csv(paste0('data/block_heritabilities_trait_', '50_irnt', '.csv'))
df$min_burden_p_no_multi[df$min_burden_p_no_multi == Inf] = NA
burden_cutoff=-log10(0.05/18524)
gwas_cutoff=-log10(5e-8)

p1 <- ggplot(df, aes(x=-log10(min_gwas_p), y=-log10(min_burden_p_no_multi))) +
  geom_point(size=1.1,color = "black",alpha=0.65) +
  geom_hline(yintercept = burden_cutoff, color = "red", linetype = "dashed") +
  geom_vline(xintercept = gwas_cutoff, linetype = "dashed",color = "red") +
  scale_x_continuous(breaks = seq(0, 352, by = 50)) +
  xlab(expression(paste(-log[10], " min p (GWAS) in LD block"))) +
  ylab(expression(paste(-log[10], " min p (burden) in LD block"))) +
  coord_cartesian(clip = "off") +
  theme_cowplot(9)

outfile="paper_figures/figS5.pdf"
save_plot(outfile, p1, base_width = 3, base_height = 3)


# Scatters for a bunch of traits
to_plot <- list()
for(idx in 1:10){
  trait_temp <- c("30050_irnt",
                  "30100_irnt",
                  "30010_irnt",
                  "30760_irnt",
                  "23105_irnt",
                  "30880_irnt",
                  "30770_irnt",
                  "30830_irnt",
                  "23106_irnt",
                  "30780_irnt")[idx]
  trait_name <- c("MCH",
                  "MPV",
                  "RBC count",
                  "HDL",
                  "BMR",
                  "Urate",
                  "IGF-1",
                  "SHBG",
                  "WBI",
                  "LDL")[idx]
  
  df <- read_csv(paste0('data/block_heritabilities_trait_', trait_temp, '.csv'))
  df$min_burden_p_no_multi[df$min_burden_p_no_multi == Inf] = NA
  p1 <- ggplot(df, aes(x=-log10(min_gwas_p), y=-log10(min_burden_p_no_multi))) +
    geom_point(size=1.1,color = "black",alpha=0.65) +
    geom_hline(yintercept = burden_cutoff, color = "red", linetype = "dashed") +
    geom_vline(xintercept = gwas_cutoff, linetype = "dashed",color = "red") +
    scale_x_continuous(breaks = seq(0, 352, by = 50)) +
    xlab(expression(paste(-log[10], " min p (GWAS) in LD block"))) +
    ylab(expression(paste(-log[10], " min p (burden) in LD block"))) +
    coord_cartesian(clip = "off") +
    theme_cowplot(9) +
    ggtitle(trait_name)
  
  to_plot[[idx]] <- p1
}

to_plot_full <- plot_grid(to_plot[[1]], to_plot[[2]], to_plot[[3]], to_plot[[4]], to_plot[[5]],
                          to_plot[[6]], to_plot[[7]], to_plot[[8]], to_plot[[9]],
                          nrow=3, ncol=3, byrow=TRUE)

save_plot('paper_figures/figS6.pdf', to_plot_full, base_width = 3*3, base_height = 3*3)




# Spearman's and bar plots
d_all_gwas_hits=fread("data/all_gwas_hits.txt")

spearmans <- list()

d_trait_names=fread("data/trait_lookup.tsv",header = F)
colnames(d_trait_names)=c("pheno","trait")

coloc <- data.frame(pheno=character(), burden_blocks=integer(), gwas_blocks=integer(), trait=character())
for (trait_temp in unique(d_all_gwas_hits$pheno)){
  df <- read_csv(paste0('data/block_heritabilities_trait_', trait_temp, '.csv'))
  trait_name <- (d_trait_names %>% filter(pheno==trait_temp))$trait
  df$min_burden_p_no_multi[df$min_burden_p_no_multi == Inf] = NA
  spearmans[[trait_temp]] <- cor(df$min_gwas_p, df$min_burden_p_no_multi, use='complete.obs', method='spearman')
  
  df <- na.omit(df, cols=c(min_gwas_p, min_burden_p_no_multi))
  sig <- df$min_burden_p_no_multi < 2.699201e-06
  num_hits <- sum(sig)
  gwas_t <- sort(df$min_gwas_p)[num_hits]
  sig_gwas <- df$min_gwas_p <= gwas_t
  overlap <- sum(sig & sig_gwas)
  
  this <- data.frame(pheno=trait_temp, burden_blocks=num_hits, gwas_blocks=overlap, trait=trait_name)
  coloc <- rbind(coloc, this)
  
}

spearman_hist <- ggplot(data.frame(unlist(spearmans)), aes(x=unlist.spearmans.)) +
  geom_histogram(bins=60) +
  xlab(expression(paste("Spearman's correlation between ", -log[10], ' min p (GWAS) and ',-log[10], ' min p (burden) in LD block'))) + 
  ylab("Number of traits") +
  theme_cowplot(9)

outfile=paste0("paper_figures/figS7.pdf")
save_plot(outfile, spearman_hist, base_width = 8, base_height = 3)



coloc$name = paste(coloc$pheno,  coloc$trait, sep=': ')
coloc$name=factor(coloc$name,levels = coloc$name)
coloc$type <- factor(coloc$name, levels = c("Burden", "Coloc"))

p1_supp<-ggplot(coloc, aes(x = name, y = burden_blocks, fill = "Burden only")) +
  geom_bar(stat = "identity", color = "black", width = 0.7) +
  geom_bar(aes(y = gwas_blocks, fill = "Burden and GWAS"), stat = "identity", color = "black", width = 0.7) +
  scale_fill_manual(values = c("Burden and GWAS" = "#1f77b4", "Burden only" = "#aec7e8"), breaks = c("", "")) +
  # scale_fill_manual(values = c("Burden and GWAS" = "#1f77b4", "Burden only" = "#aec7e8"),
  #                   breaks = c("Burden only", "Burden and GWAS")) +
  labs(x = "", y = "Number of LD blocks") +
  theme_cowplot(7) +
  coord_flip() + 
  theme(# axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = c(0.3, 0.9),
    legend.title = element_blank())

outfile=paste0("paper_figures/figS4.pdf")
save_plot(outfile, p1_supp, base_width = 2.25*4, base_height = 2*5)









#########################
### Figure 3B-D, 6B ###
#########################


d_indep_traits=fread("data/indep_traits.txt")

d_s_het=fread("data/zeng_et_al_v3.tsv")
d_gamma_unbiased=fread("data/unbiased_gamma_sq_by_trait_MAF1_hgncIDs.txt")
d_gamma_hats=fread("data/gamma_hats_by_trait_MAF1_hgncIDs.txt")
d_z_scores=fread("data/zscores_by_trait_MAF1_hgncIDs.txt")
d_SEs=fread("data/SE_hats_by_trait_MAF1_hgncIDs.txt")
d_plof = read.csv('data/plof_MAF1_hgncIDs.csv')

###

set.seed(0)

traits=sort(unique(d_indep_traits$trait)) 
traits=traits[!(traits %in% c("20127_irnt","2139_irnt"))]

s_het_bin_size = 100
colnames(d_s_het)[2]="hgnc_id"

d_s_het=d_s_het[!is.na(d_s_het$post_mean),]
d_s_het=d_s_het[d_s_het$hgnc_id %in% d_gamma_hats$hgnc_id,]
d_s_het$s_het_cat=ntile(d_s_het$post_mean,s_het_bin_size)


#

d_gamma2_matrix= d_gamma_unbiased %>% dplyr::select(-hgnc_id)  %>% dplyr::select(traits)
d_SE_matrix= d_SEs %>% dplyr::select(-hgnc_id) %>% dplyr::select(traits)
d_z_matrix= d_z_scores %>% dplyr::select(-hgnc_id) %>% dplyr::select(traits)
d_z_matrix=replace(d_z_matrix, is.na(d_z_matrix), 1) #replace NAs with 1

d_test_gamma <- ifelse(abs(d_gamma2_matrix) < 4, 1, 1)
d_test_se <- ifelse(abs(d_SE_matrix) < 4, 1, 1)
d_test_z <- ifelse(abs(d_z_matrix) < 4, 1, 1)

d_data=(d_z_scores %>% dplyr::select(hgnc_id))
d_data$mean_z2=rowMeans(d_z_matrix^2,na.rm = T)
d_data$mean_gamma2=rowMeans(d_gamma2_matrix,na.rm = T)
d_data$mean_SE2=rowMeans(d_SE_matrix^2,na.rm = T)
d_data$count_traits_gamma=rowSums(d_test_gamma,na.rm = T)
d_data$count_traits_se=rowSums(d_test_se,na.rm = T)
d_data$count_traits_z=rowSums(d_test_z,na.rm = T)

d_data=left_join(d_data, d_s_het %>% dplyr::select(hgnc_id,exp_lof,post_mean,s_het_cat))
d_data=left_join(d_data, d_plof)
d_data=d_data[!is.na(d_data$post_mean),]

data_sum <- d_data %>%
  group_by(s_het_cat) %>%
  summarize(mean_gamma_sq = mean(mean_gamma2, na.rm = TRUE),
            mean_se_sq = mean(mean_SE2, na.rm = TRUE),
            mean_z_sq=mean(mean_z2, na.rm = TRUE),
            mean_s = mean(post_mean, na.rm = TRUE),
            sd_gamma_sq = sd(mean_gamma2, na.rm = TRUE),
            sd_se_sq = sd(mean_SE2, na.rm = TRUE),
            sd_z_sq = sd(mean_z2, na.rm = TRUE),
            count_gamma_sq = sum(count_traits_gamma),
            count_se_sq = sum(count_traits_se),
            count_z_sq = sum(count_traits_z),
            mean_plof = mean(effect_allele_frequency, na.rm=TRUE)
  )

data_sum$se_gamma_sq=data_sum$sd_gamma_sq/sqrt(data_sum$count_gamma_sq)
data_sum$se_se_sq=data_sum$sd_se_sq/sqrt(data_sum$count_se_sq)
data_sum$se_z_sq=data_sum$sd_z_sq/sqrt(data_sum$count_z_sq)

data_enrich <- d_z_matrix^2 - 1
data_enrich$hgnc_id <- (d_z_scores %>% dplyr::select(hgnc_id))
data_enrich <- left_join(data_enrich, d_s_het %>% dplyr::select(hgnc_id,exp_lof,post_mean,s_het_cat))
data_enrich <- data_enrich[!is.na(data_enrich$post_mean),]
traits <- colnames(data_enrich)[grepl('_irnt', colnames(data_enrich))]
data_enrich <- data_enrich %>% mutate_at(traits, function(x) x / mean(x))
data_enrich_var <- data_enrich %>% group_by(s_het_cat) %>% summarize_if(function(x) is.numeric(x) , var)
data_enrich <- data_enrich %>% group_by(s_het_cat) %>% summarize_if(function(x) is.numeric(x) , mean)
data_enrich$avg_h2_enrich <- rowMeans(data_enrich[traits] / data_enrich_var[traits]) / rowMeans(1/data_enrich_var[traits])
data_enrich$avg_h2_enrich <- data_enrich$avg_h2_enrich / mean(data_enrich$avg_h2_enrich)

#########
span_temp=0.8

line_alpha=0.8
point_alpha=0.2

#########

color_temp="purple"
p0 <- ggplot(data_enrich, aes(x=post_mean, y=avg_h2_enrich)) +
  geom_point(size=2,alpha=point_alpha,color=color_temp) + 
  geom_smooth(method = "loess", span=span_temp,se = F, color = alpha(color_temp, line_alpha)) +
  xlab(substitute(paste("Mean s"[het], " in bin"))) +
  ylab(expression("Heritability enrichment")) +
  scale_x_continuous(trans='log', breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1)) +
  theme_cowplot(9)

outfile=paste0("paper_figures/fig6d.pdf")
save_plot(outfile, p0,base_width = 2.5, base_height = 2)

#
color_temp="#4169E1"
p1 <- ggplot(data_sum, aes(x=mean_s, y=mean_gamma_sq)) +
  geom_point(size=2,alpha=point_alpha,color=color_temp) + 
  geom_smooth(method = "loess", span=span_temp,se = F, color = alpha(color_temp, line_alpha)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  xlab(substitute(paste("Mean s"[het], " in bin"))) +
  ylab(expression(paste("Mean(" * gamma^2 * ")" ~ "across traits"))) +
  theme_cowplot(9)

outfile=paste0("paper_figures/fig3c.pdf")
save_plot(outfile, p1,base_width = 2.5, base_height = 2)

#

color_temp="#DC143C"
p2 <- ggplot(data_sum, aes(x=mean_s, y=mean_plof)) +
  geom_point(size=2,alpha=point_alpha,color=color_temp) + 
  geom_smooth(method = "loess", span=span_temp,se = F, color = alpha(color_temp, line_alpha)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  xlab(substitute(paste("Mean s"[het], " in bin"))) +
  ylab(expression(paste("Mean p"[LoF], " in bin"))) +
  theme_cowplot(9)

outfile=paste0("paper_figures/fig3b.pdf")
save_plot(outfile, p2,base_width = 2.5, base_height = 2)

#

color_range <- c(-3.4,0)
color_breaks <- c(-3, -2, -1, 0)
color_labels <- expression(10^-3,10^-2,10^-1,10^0)

color_scheme <- scales::col_numeric(
  palette = rev(c("black", "red", "orange", "yellow")),
  domain = color_range
)

color_temp="#333333"

p3 <- ggplot(data_sum, aes(x=mean_gamma_sq, y=mean_z_sq,colour = log10(mean_s))) +
  geom_point(size=2,alpha=0.8) + 
  geom_smooth(method = "loess", span=span_temp,se = F, color = alpha(color_temp, line_alpha)) +
  xlab(expression(paste("Mean(" * gamma^2 * ")" ~ "across traits"))) +
  ylab(expression("Mean(" * z^2 * ")" ~ "across traits")) +
  scale_color_gradient(name = substitute(paste("Mean s"[het])), limits = color_range, breaks = color_breaks, labels = color_labels,
                       low = color_scheme(-3.4), high = color_scheme(0),
                       guide = guide_colorbar(label.theme = element_text(hjust = 0,size=9)))+  
  scale_x_continuous()+
  theme_cowplot(9)

outfile=paste0("paper_figures/fig3d.pdf")
save_plot(outfile, p3,base_width = 3, base_height = 2)







#########################
### Figure S8 ###
#########################

d_indep_traits=fread("data/indep_traits.txt")

d_s_het=fread("data/zeng_et_al_v3.tsv")
d_gamma_unbiased=fread("data/unbiased_gamma_sq_by_trait_MAF1_hgncIDs.txt")
d_gamma_hats=fread("data/gamma_hats_by_trait_MAF1_hgncIDs.txt")
d_z_scores=fread("data/zscores_by_trait_MAF1_hgncIDs.txt")
d_SEs=fread("data/SE_hats_by_trait_MAF1_hgncIDs.txt")
d_plof = read.csv('data/plof_MAF1_hgncIDs.csv')

###

set.seed(0)

traits=sort(unique(d_indep_traits$trait)) 
traits=traits[!(traits %in% c("20127_irnt","2139_irnt"))]

s_het_bin_size = 100
colnames(d_s_het)[2]="hgnc_id"

d_s_het=d_s_het[!is.na(d_s_het$prior_mean),]
d_s_het=d_s_het[d_s_het$hgnc_id %in% d_gamma_hats$hgnc_id,]
d_s_het$s_het_cat=ntile(d_s_het$prior_mean,s_het_bin_size)


#

d_gamma2_matrix= d_gamma_unbiased %>% dplyr::select(-hgnc_id)  %>% dplyr::select(traits)
d_SE_matrix= d_SEs %>% dplyr::select(-hgnc_id) %>% dplyr::select(traits)
d_z_matrix= d_z_scores %>% dplyr::select(-hgnc_id) %>% dplyr::select(traits)
d_z_matrix=replace(d_z_matrix, is.na(d_z_matrix), 1) #replace NAs with 1

d_test_gamma <- ifelse(abs(d_gamma2_matrix) < 4, 1, 1)
d_test_se <- ifelse(abs(d_SE_matrix) < 4, 1, 1)
d_test_z <- ifelse(abs(d_z_matrix) < 4, 1, 1)

d_data=(d_z_scores %>% dplyr::select(hgnc_id))
d_data$mean_z2=rowMeans(d_z_matrix^2,na.rm = T)
d_data$mean_gamma2=rowMeans(d_gamma2_matrix,na.rm = T)
d_data$mean_SE2=rowMeans(d_SE_matrix^2,na.rm = T)
d_data$count_traits_gamma=rowSums(d_test_gamma,na.rm = T)
d_data$count_traits_se=rowSums(d_test_se,na.rm = T)
d_data$count_traits_z=rowSums(d_test_z,na.rm = T)

d_data=left_join(d_data, d_s_het %>% dplyr::select(hgnc_id,exp_lof,prior_mean,s_het_cat))
d_data=left_join(d_data, d_plof)
d_data=d_data[!is.na(d_data$prior_mean),]

data_sum <- d_data %>%
  group_by(s_het_cat) %>%
  summarize(mean_gamma_sq = mean(mean_gamma2, na.rm = TRUE),
            mean_se_sq = mean(mean_SE2, na.rm = TRUE),
            mean_z_sq=mean(mean_z2, na.rm = TRUE),
            mean_s = mean(prior_mean, na.rm = TRUE),
            sd_gamma_sq = sd(mean_gamma2, na.rm = TRUE),
            sd_se_sq = sd(mean_SE2, na.rm = TRUE),
            sd_z_sq = sd(mean_z2, na.rm = TRUE),
            count_gamma_sq = sum(count_traits_gamma),
            count_se_sq = sum(count_traits_se),
            count_z_sq = sum(count_traits_z),
            mean_plof = mean(effect_allele_frequency, na.rm=TRUE)
  )

data_sum$se_gamma_sq=data_sum$sd_gamma_sq/sqrt(data_sum$count_gamma_sq)
data_sum$se_se_sq=data_sum$sd_se_sq/sqrt(data_sum$count_se_sq)
data_sum$se_z_sq=data_sum$sd_z_sq/sqrt(data_sum$count_z_sq)

data_enrich <- d_z_matrix^2 - 1
data_enrich$hgnc_id <- (d_z_scores %>% dplyr::select(hgnc_id))
data_enrich <- left_join(data_enrich, d_s_het %>% dplyr::select(hgnc_id,exp_lof,prior_mean,s_het_cat))
data_enrich <- data_enrich[!is.na(data_enrich$prior_mean),]
traits <- colnames(data_enrich)[grepl('_irnt', colnames(data_enrich))]
data_enrich <- data_enrich %>% mutate_at(traits, function(x) x / mean(x))
data_enrich_var <- data_enrich %>% group_by(s_het_cat) %>% summarize_if(function(x) is.numeric(x) , var)
data_enrich <- data_enrich %>% group_by(s_het_cat) %>% summarize_if(function(x) is.numeric(x) , mean)
data_enrich$avg_h2_enrich <- rowMeans(data_enrich[traits] / data_enrich_var[traits]) / rowMeans(1/data_enrich_var[traits])
data_enrich$avg_h2_enrich <- data_enrich$avg_h2_enrich / mean(data_enrich$avg_h2_enrich)

#########
span_temp=0.8

line_alpha=0.8
point_alpha=0.2


color_temp="#DC143C"
p2 <- ggplot(data_sum, aes(x=mean_s, y=mean_plof)) +
  geom_point(size=2,alpha=point_alpha,color=color_temp) + 
  geom_smooth(method = "loess", span=span_temp,se = F, color = alpha(color_temp, line_alpha)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  xlab(substitute(paste("Mean prior mean s"[het], " in bin"))) +
  ylab(expression(paste("Mean p"[LoF], " in bin"))) +
  theme_cowplot(9)

outfile=paste0("paper_figures/figS8.pdf")
save_plot(outfile, p2,base_width = 4.5, base_height = 4)









#########################
### Figure 3E ###
#########################

squish_t <- 30
squish_factor <- 0.75
squish_trans <- function() {
  trans_new(
    name = "squish",
    transform = function(x) ifelse(x <= squish_t, x, log(x-squish_t + 1)/squish_factor + squish_t),
    inverse = function(x) ifelse(x <= squish_t, x, exp((x-squish_t)*squish_factor) - 1 + squish_t)
  )
}



convert_name_to_hgnc <- function(gene_ids,d_conv){
  
  dx=d_conv[d_conv$Symbol %in% gene_ids,]
  dx=dx[!duplicated(dx$Symbol),]
  
  dx=dx[!duplicated(dx$hgnc_id),]
  
  colnames(dx)=c("gene","hgnc_id")
  return(dx)
  
}

d_z_scores <- fread("data/zscores_by_trait_MAF1_hgncIDs.txt")

# 30040_irnt  Erythroidcells
# 30240_irnt  Erythroidcells
# 30210_irnt  T-cells
# 30120_irnt  T-cells
# 50_irnt bone
# 3148_irnt   bone
# 30740_irnt  pancreas
# 30620_irnt  liver
# 30700_irnt  liver

d_full <- data.frame(gene=character(), hgnc_id=character(), num_tissues=integer(), z=double(), expr_group=character(), p_val=double())
percs <- c(0, 0, 0, 0)
plots <- list()
for(idx in 1:9){
  t <- c('30040_irnt', '30240_irnt', '30210_irnt',
         '30120_irnt', '50_irnt', '3148_irnt',
         '30740_irnt', '30620_irnt', '30700_irnt')[idx]
  cell_type <- c('Erythroidcells', 'Erythroidcells', 'T-cells',
                 'T-cells', 'bone', 'bone',
                 'pancreas', 'liver', 'liver')[idx]
  d_conv <- fread('data/all_genes.hgnc_id.txt')
  exp_name <- paste('data/', cell_type, '_TPM10Genes_expression_specificity.txt', sep='')
  d_exp <- read.table(exp_name, col.names=c('gene', 'num_tissues'))
  d_conv <- convert_name_to_hgnc(d_exp$gene, d_conv)
  d_exp <- left_join(d_conv, d_exp, by='gene')
  these_z_scores <- d_z_scores %>% select(all_of(c('hgnc_id', t)))
  colnames(these_z_scores)[2] = 'z'
  this_full <- left_join(d_exp, these_z_scores, by='hgnc_id')
  this_full$trait <- t
  this_full <- drop_na(this_full, 'z')
  this_full$num_tissues <- as.numeric(this_full$num_tissues)
  this_full$expr_group <- cut(this_full$num_tissues, breaks = c(-1,0.041, 0.061, 0.088, 0.149, 2), labels = c("group1", "group2", "group3", "group4", "group5"))
  this_full$p_val <- pchisq(this_full$z^2, 1, lower.tail=FALSE, log.p = TRUE)
  plots[[idx]] <- ggplot(this_full, aes(sample=-p_val/log(10), colour=expr_group)) +
    stat_qq(distribution = function(p) qexp(p, rate=log(10)), alpha=1, size=2) + # stat_qq_line(distribution = function(p) qexp(p, rate=1)) +
    xlab(expression("Expected "-log[10](p))) +
    ylab(expression("Observed "-log[10](p))) +
    ggtitle(paste(t, cell_type)) +
    geom_abline() +
    guides(colour="none") +
    scale_y_continuous(trans = squish_trans(), breaks=c(0, 10, 20, 30, 60, 300)) +
    theme_cowplot(9)
  d_full <- rbind(d_full, this_full)
}


p1 <- ggplot(d_full, aes(sample=-p_val/log(10), colour=expr_group)) +
  stat_qq(distribution = function(p) qexp(p, rate=log(10)), alpha=0.5, size=2) + # stat_qq_line(distribution = function(p) qexp(p, rate=1)) +
  xlab(expression("Expected "-log[10](p))) +
  ylab(expression("Observed "-log[10](p))) +
  geom_abline() + 
  guides(color=guide_legend(reverse=TRUE)) +
  scale_y_continuous(trans = squish_trans(), breaks=c(0, 10, 20, 30, 60, 300)) +
  scale_color_discrete(name='Expression\nspecificity bin', labels=c('Least specific', '', '', '', 'Most specific')) +
  theme_cowplot(9)

outfile="paper_figures/fig3e.pdf"
save_plot(outfile, p1,base_width = 4.2, base_height = 2)











#########################
### Figure S9 ###
#########################



d_data=fread("data/props_by_psi_expr_jeff.txt")
d_data <- separate(d_data, trait_cell, into = c("trait", "tissue"), sep = "_irnt_")

tissues=sort(unique(d_data$tissue))
tissue_names=c("Bone","Erythroid","Liver","Pancreas","Tcell")
d_map=data.frame(tissue=tissues,tissue_name=tissue_names)
d_data=left_join(d_data,d_map)

traits=sort(unique(d_data$trait))
#trait_names=c("Mean corpuscular volume","Lymphocyte count","Eosinophil count","Reticulocyte percentage", "Alanine aminotransferase", "Creatinine", "Glucose", "Heel bone mineral density", "Height")
trait_names=c("MCV","Lymphocyte count","Eosinophil count","Reticulocyte percentage", "ALT", "Creatinine", "Glucose", "Heel BMD", "Height")
d_map=data.frame(trait=traits,trait_name=trait_names)
d_data=left_join(d_data,d_map)

d_data$trait_cell=paste0(d_data$trait_name,"-",d_data$tissue_name)

###meta-analysis reg_coefs

d_data$weight <- 1 / (d_data$reg_coef_se)^2

iter=0
for (group_temp in c("group1", "group2", "group3", "group4", "group5")){
  iter=iter+1
  d_temp=d_data[d_data$expr_group==group_temp,]
  weighted_mean <- sum(d_temp$weight * d_temp$reg_coef) / sum(d_temp$weight)
  weighted_se <- sqrt(1 / sum(d_temp$weight))
  if (group_temp == "group1"){weighted_mean=0; weighted_se=0}
  d_sum_temp=data.frame(expr_group= group_temp,reg_coef=weighted_mean,reg_coef_se=weighted_se,trait_cell="Average")
  if (iter==1){d_sum=d_sum_temp}
  if (iter>1){d_sum=rbind(d_sum,d_sum_temp)}
}

########

d_data$expr_group=factor(d_data$expr_group, levels = (c("group1", "group2", "group3", "group4", "group5")))
d_sum$expr_group=factor(d_sum$expr_group, levels = (c("group1", "group2", "group3", "group4", "group5")))

d_all=rbind(d_sum,d_data %>% select(colnames(d_sum)))
d_all$trait_cell <- factor(d_all$trait_cell,levels = c(sort(unique(d_data$trait_cell)),"Average")  )

########

num_levels <- length(levels(d_all$trait_cell))
custom_palette <- brewer.pal(num_levels - 1, "Set3")
custom_palette <- c(custom_palette, "black")

b_point = 0.4
squish_power = 0.1
squish_const = 0.2
squish_trans <- function() {
  trans_new(
    name = "squish",
    transform = function(x) sign(x)*((abs(x/b_point)+squish_const)^(squish_power)*squish_const^(squish_power-1)/squish_power - squish_const^(2*squish_power-1)/squish_power),
    inverse = function(x) sign(x)*b_point*(((squish_power*abs(x)+squish_const^(2*squish_power-1))/squish_const^(squish_power-1))^(1/squish_power) - squish_const)
  )
}

breaks <- c(0, 0.4, 0.8, 1.2, 1.6, 2.0)
labels <- c(0, 0.4, 0.8, 1.2, 1.6, 2.0)

p0 <- ggplot() +
  geom_line(data=d_all,aes(x=expr_group, y=reg_coef,group=trait_cell,color=trait_cell),linewidth=0.5) +
  geom_point(data=d_sum,aes(x=expr_group, y=reg_coef),size=2.2) +
  geom_errorbar(data=d_sum,aes(x=expr_group, y=reg_coef,ymin=reg_coef-2*reg_coef_se, ymax=reg_coef+2*reg_coef_se), width=0, position=position_dodge(0.0))+
  geom_line(data=d_sum,aes(x=expr_group, y=reg_coef,group = 1),linewidth=1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black",linewidth=0.5) +
  xlab("Expression specificity bin") +
  ylab(expression(paste("Effect on ", z^2))) +
  scale_x_discrete(labels = c("Bin 1", "Bin 2", "Bin 3","Bin 4","Bin 5")) +
  labs(color = "Trait-Tissue pair") +
  scale_color_manual(name = "Trait-Tissue pair", values = custom_palette) +
  scale_y_continuous(trans = squish_trans(), breaks = breaks, limits=c(-0.1, 2.2)) +
  theme_cowplot(9)+
  guides(color = "none")

outfile="paper_figures/figS9b.pdf"
save_plot(outfile, p0,base_width = 3.5, base_height = 3)

####

d_data=fread("data/props_by_psi_expr_jeff_no_correct.txt")
d_data <- separate(d_data, trait_cell, into = c("trait", "tissue"), sep = "_irnt_")

tissues=sort(unique(d_data$tissue))
tissue_names=c("Bone","Erythroid","Liver","Pancreas","Tcell")
d_map=data.frame(tissue=tissues,tissue_name=tissue_names)
d_data=left_join(d_data,d_map)

traits=sort(unique(d_data$trait))
#trait_names=c("Mean corpuscular volume","Lymphocyte count","Eosinophil count","Reticulocyte percentage", "Alanine aminotransferase", "Creatinine", "Glucose", "Heel bone mineral density", "Height")
trait_names=c("MCV","Lymphocyte count","Eosinophil count","Reticulocyte percentage", "ALT", "Creatinine", "Glucose", "Heel BMD", "Height")
d_map=data.frame(trait=traits,trait_name=trait_names)
d_data=left_join(d_data,d_map)

d_data$trait_cell=paste0(d_data$trait_name,"-",d_data$tissue_name)

###meta-analysis reg_coefs

d_data$weight <- 1 / (d_data$reg_coef_se)^2

iter=0
for (group_temp in c("group1", "group2", "group3", "group4", "group5")){
  iter=iter+1
  d_temp=d_data[d_data$expr_group==group_temp,]
  weighted_mean <- sum(d_temp$weight * d_temp$reg_coef) / sum(d_temp$weight)
  weighted_se <- sqrt(1 / sum(d_temp$weight))
  if (group_temp == "group1"){weighted_mean=0; weighted_se=0}
  d_sum_temp=data.frame(expr_group= group_temp,reg_coef=weighted_mean,reg_coef_se=weighted_se,trait_cell="Average")
  if (iter==1){d_sum=d_sum_temp}
  if (iter>1){d_sum=rbind(d_sum,d_sum_temp)}
}

########

d_data$expr_group=factor(d_data$expr_group, levels = (c("group1", "group2", "group3", "group4", "group5")))
d_sum$expr_group=factor(d_sum$expr_group, levels = (c("group1", "group2", "group3", "group4", "group5")))

d_all=rbind(d_sum,d_data %>% select(colnames(d_sum)))
d_all$trait_cell <- factor(d_all$trait_cell,levels = c(sort(unique(d_data$trait_cell)),"Average")  )

########

num_levels <- length(levels(d_all$trait_cell))
custom_palette <- brewer.pal(num_levels - 1, "Set3")
custom_palette <- c(custom_palette, "black")

b_point = 0.4
squish_power = 0.1
squish_const = 0.2
squish_trans <- function() {
  trans_new(
    name = "squish",
    transform = function(x) sign(x)*((abs(x/b_point)+squish_const)^(squish_power)*squish_const^(squish_power-1)/squish_power - squish_const^(2*squish_power-1)/squish_power),
    inverse = function(x) sign(x)*b_point*(((squish_power*abs(x)+squish_const^(2*squish_power-1))/squish_const^(squish_power-1))^(1/squish_power) - squish_const)
  )
}

breaks <- c(0, 0.4, 0.8, 1.2, 1.6, 2.0)
labels <- c(0, 0.4, 0.8, 1.2, 1.6, 2.0)

p0 <- ggplot() +
  geom_line(data=d_all,aes(x=expr_group, y=reg_coef,group=trait_cell,color=trait_cell),linewidth=0.5) +
  geom_point(data=d_sum,aes(x=expr_group, y=reg_coef),size=2.2) +
  geom_errorbar(data=d_sum,aes(x=expr_group, y=reg_coef,ymin=reg_coef-2*reg_coef_se, ymax=reg_coef+2*reg_coef_se), width=0, position=position_dodge(0.0))+
  geom_line(data=d_sum,aes(x=expr_group, y=reg_coef,group = 1),linewidth=1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black",linewidth=0.5) +
  xlab("Expression specificity bin") +
  ylab(expression(paste("Effect on ", z^2))) +
  scale_x_discrete(labels = c("Bin 1", "Bin 2", "Bin 3","Bin 4","Bin 5")) +
  labs(color = "Trait-Tissue pair") +
  scale_color_manual(name = "Trait-Tissue pair", values = custom_palette) +
  scale_y_continuous(trans = squish_trans(), breaks = breaks, limits=c(-0.1, 2.2)) +
  theme_cowplot(9)+
  guides(color = "none")



outfile="paper_figures/figS9a.pdf"
save_plot(outfile, p0,base_width = 3.5, base_height = 3)







#########################
### Figure 4A ###
#########################

coding_df <- data.table()
coding_df$x <- c(0., 1.)
coding_df$y <- c(1., 1.)
coding_df$psi_v <- c(NA, NA)


psi_curve_df <- data.table(psi_g=vector(), context_factor=vector(), psi_v=vector())
for(this_psi_v in seq(0.01, 1.0, length.out=20)){
  xs <- exp(seq(log(this_psi_v/10), log(1), length.out=30))
  ys <- this_psi_v/xs
  this_df <- data.table(psi_g=xs, context_factor=ys, psi_v=rep(this_psi_v, length(xs)))
  psi_curve_df = rbind(psi_curve_df, this_df)
}

p1 <- ggplot(psi_curve_df, aes(x=psi_g, y=context_factor, group=psi_v)) + 
  geom_line(aes(color=psi_v)) +
  scale_color_viridis(direction = -1) +
  xlab(expression(Psi["G"] ~ ", trait specificity of the gene")) +
  ylab("Trait specificity of variant relative to gene") +
  geom_line(data=data.frame(coding_df), aes(x=x, y=y), linetype=2, linewidth=1) +
  guides(color = "none")

outfile="paper_figures/fig4a.pdf"
save_plot(outfile, p1, base_width = 3, base_height = 3)











#########################
### Figure 4B and S10 ###
#########################

d_ldsc = fread('data/all_traits_spcificityBIN_tau_jeff_bins.txt')

d_ldsc <- separate(d_ldsc, annot, into = c("tissue", "covar", "group"), sep = "_")
d_data=d_ldsc[d_ldsc$covar=="specificity",]

tissues=sort(unique(d_data$tissue))
tissue_names=c("Bone","Erythroid","Liver","Pancreas","Tcell")
d_map=data.frame(tissue=tissues,tissue_name=tissue_names)
d_data=left_join(d_data,d_map)

traits=sort(unique(d_data$trait))
trait_names=c("MCV","Lymphocyte count","Eosinophil count","Reticulocyte percentage", "ALT", "Creatinine", "Glucose", "Heel BMD", "Height")
d_map=data.frame(trait=traits,trait_name=trait_names)
d_data=left_join(d_data,d_map)

d_data$trait_cell=paste0(d_data$trait_name,"-",d_data$tissue_name)

###meta-analysis

d_data$weight <- 1 / (d_data$tau_se)^2

iter=0
for (group_temp in c("GROUP1", "GROUP2", "GROUP3", "GROUP4", "GROUP5")){
  iter=iter+1
  d_temp=d_data[d_data$group==group_temp,]
  weighted_mean <- sum(d_temp$weight * d_temp$tau) / sum(d_temp$weight)
  weighted_se <- sqrt(1 / sum(d_temp$weight))
  d_sum_temp=data.frame(group= group_temp,tau=weighted_mean,tau_se=weighted_se,trait_cell="Average")
  if (iter==1){d_sum=d_sum_temp}
  if (iter>1){d_sum=rbind(d_sum,d_sum_temp)}
}

########

d_data$group=factor(d_data$group, levels = c("GROUP1", "GROUP2", "GROUP3", "GROUP4", "GROUP5"))
d_sum$group=factor(d_sum$group, levels = c("GROUP1", "GROUP2", "GROUP3", "GROUP4", "GROUP5"))

d_all=rbind(d_sum,d_data %>% select(colnames(d_sum)))
d_all$trait_cell <- factor(d_all$trait_cell,levels = c(sort(unique(d_data$trait_cell)),"Average")  )

########

num_levels <- length(levels(d_all$trait_cell))
custom_palette <- brewer.pal(num_levels - 1, "Set3")
custom_palette <- c(custom_palette, "black")

###


b_point = 2.5
squish_power = 0.1
squish_const = 0.1
squish_trans <- function() {
  trans_new(
    name = "squish",
    transform = function(x) sign(x)*((abs(x/b_point)+squish_const)^(squish_power)*squish_const^(squish_power-1)/squish_power - squish_const^(2*squish_power-1)/squish_power),
    inverse = function(x) sign(x)*b_point*(((squish_power*abs(x)+squish_const^(2*squish_power-1))/squish_const^(squish_power-1))^(1/squish_power) - squish_const)
  )
}


breaks <- c(-2, -1, 0, 1, 2, 3, 4)
labels <- c(-2, -1, 0, 1, 2, 3, 4)

###

p0 <- ggplot() +
  geom_line(data=d_all,aes(x=group, y=tau,group=trait_cell,color=trait_cell),linewidth=0.5) +
  geom_point(data=d_sum,aes(x=group, y=tau),size=2.2) +
  geom_errorbar(data=d_sum,aes(x=group, y=tau,ymin=tau-2*tau_se, ymax=tau+2*tau_se), width=0, position=position_dodge(0.0))+
  geom_line(data=d_sum,aes(x=group, y=tau,group = 1),linewidth=1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black",linewidth=0.5) +
  xlab("Expression specificity bin") +
  ylab(expression("LDSC" ~ tau ~ ", coding variants")) +
  scale_x_discrete(labels = c("Bin 1", "Bin 2", "Bin 3", "Bin 4", "Bin 5")) +
  labs(color = "Trait-Tissue pair") +
  scale_color_manual(name = "Trait-Tissue pair", values = custom_palette) +
  theme_cowplot(9) +
  guides(color = "none") +
  scale_y_continuous(trans = squish_trans(), breaks = breaks, labels = labels)


outfile="paper_figures/fig4b.pdf"
save_plot(outfile, p0,base_width = 2.3, base_height = 2.1)

d_ldsc = fread('data/all_traits_intensityBIN_tau_jeff_bins.txt')

d_ldsc <- separate(d_ldsc, annot, into = c("tissue", "covar", "group"), sep = "_")
d_data=d_ldsc[d_ldsc$covar=="TPM",]

tissues=sort(unique(d_data$tissue))
tissue_names=c("Bone","Erythroid","Liver","Pancreas","Tcell")
d_map=data.frame(tissue=tissues,tissue_name=tissue_names)
d_data=left_join(d_data,d_map)

traits=sort(unique(d_data$trait))
trait_names=c("MCV","Lymphocyte count","Eosinophil count","Reticulocyte percentage", "ALT", "Creatinine", "Glucose", "Heel BMD", "Height")
d_map=data.frame(trait=traits,trait_name=trait_names)
d_data=left_join(d_data,d_map)

d_data$trait_cell=paste0(d_data$trait_name,"-",d_data$tissue_name)

###meta-analysis

d_data$weight <- 1 / (d_data$tau_se)^2

iter=0
for (group_temp in c("GROUP1", "GROUP2", "GROUP3", "GROUP4")){
  iter=iter+1
  d_temp=d_data[d_data$group==group_temp,]
  weighted_mean <- sum(d_temp$weight * d_temp$tau) / sum(d_temp$weight)
  weighted_se <- sqrt(1 / sum(d_temp$weight))
  d_sum_temp=data.frame(group= group_temp,tau=weighted_mean,tau_se=weighted_se,trait_cell="Average")
  if (iter==1){d_sum=d_sum_temp}
  if (iter>1){d_sum=rbind(d_sum,d_sum_temp)}
}

########

d_data$group=factor(d_data$group, levels = rev(c("GROUP1", "GROUP2", "GROUP3", "GROUP4")))
d_sum$group=factor(d_sum$group, levels = rev(c("GROUP1", "GROUP2", "GROUP3", "GROUP4")))

d_all=rbind(d_sum,d_data %>% select(colnames(d_sum)))
d_all$trait_cell <- factor(d_all$trait_cell,levels = c(sort(unique(d_data$trait_cell)),"Average")  )

p0 <- ggplot() +
  geom_line(data=d_all,aes(x=group, y=tau,group=trait_cell,color=trait_cell),linewidth=0.5) +
  geom_point(data=d_sum,aes(x=group, y=tau),size=2.2) +
  geom_errorbar(data=d_sum,aes(x=group, y=tau,ymin=tau-2*tau_se, ymax=tau+2*tau_se), width=0, position=position_dodge(0.0))+
  geom_line(data=d_sum,aes(x=group, y=tau,group = 1),linewidth=1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black",linewidth=0.5) +
  xlab("Expression level bin") +
  ylab(expression("LDSC" ~ tau ~ ", coding variants")) +
  scale_x_discrete(labels = c("Bin 2", "Bin 3", "Bin 4", "Bin 5")) +
  labs(color = "Trait-Tissue pair") +
  scale_color_manual(name = "Trait-Tissue pair", values = custom_palette) +
  theme_cowplot(9) +
  guides(color = "none") +
  scale_y_continuous(trans = squish_trans(), breaks = breaks, labels = labels)


outfile="paper_figures/figS10.pdf"
save_plot(outfile, p0,base_width = 3.5, base_height = 3)

















#########################
### Figure 4C and S11 ###
#########################


d_ldsc=fread("data/ATAC_specificity_LDSC_Mineto.txt")

d_ldsc <- separate(d_ldsc, annot, into = c("tissue", "x", "group"), sep = "_")
d_ldsc <- separate(d_ldsc, x, into = c("covar", "y"), sep = "-") %>% select(-y)

d_data=d_ldsc[d_ldsc$covar=="specificity",]

tissues=sort(unique(d_data$tissue))
traits=sort(unique(d_data$trait))
trait_names=c("MCV","Lymphocyte count","Eosinophil count","Reticulocyte %", "ALT", "Creatinine", "Glucose", "Heel BMD", "Height")
d_map=data.frame(trait=traits,trait_name=trait_names)

d_data=left_join(d_data,d_map)
d_data$trait_cell=paste0(d_data$trait_name,"-",d_data$tissue)

###meta-analysis

d_data$weight <- 1 / (d_data$tau_se)^2

iter=0
for (group_temp in c("GROUP1", "GROUP2", "GROUP3", "GROUP4", "GROUP5")){
  iter=iter+1
  d_temp=d_data[d_data$group==group_temp,]
  weighted_mean <- sum(d_temp$weight * d_temp$tau) / sum(d_temp$weight)
  weighted_se <- sqrt(1 / sum(d_temp$weight))
  d_sum_temp=data.frame(group= group_temp,tau=weighted_mean,tau_se=weighted_se,trait_cell="Average")
  if (iter==1){d_sum=d_sum_temp}
  if (iter>1){d_sum=rbind(d_sum,d_sum_temp)}
}

########

d_data$group=factor(d_data$group, levels = rev(c("GROUP1", "GROUP2", "GROUP3", "GROUP4", "GROUP5")))
d_sum$group=factor(d_sum$group, levels = rev(c("GROUP1", "GROUP2", "GROUP3", "GROUP4", "GROUP5")))

d_all=rbind(d_sum,d_data %>% select(colnames(d_sum)))
d_all$trait_cell <- factor(d_all$trait_cell,levels = c(sort(unique(d_data$trait_cell)),"Average")  )

########

b_point = 0.65
squish_power = 0.1
squish_const = 0.2
squish_trans <- function() {
  trans_new(
    name = "squish",
    transform = function(x) sign(x)*((abs(x)+squish_const)^(squish_power)*squish_const^(squish_power-1)/squish_power - squish_const^(2*squish_power-1)/squish_power),
    inverse = function(x) sign(x)*(((squish_power*abs(x)+squish_const^(2*squish_power-1))/squish_const^(squish_power-1))^(1/squish_power) - squish_const)
  )
}

num_levels <- length(levels(d_all$trait_cell))
custom_palette <- brewer.pal(num_levels - 1, "Set3")
custom_palette <- c(custom_palette, "black")


breaks <- c(-4, -3, -2, -1, 0, 1, 2)
labels <- c(-4, -3, -2, -1, 0, 1, 2)

###

p0 <- ggplot() +
  geom_line(data=d_all,aes(x=group, y=tau,group=trait_cell,color=trait_cell),linewidth=0.5) +
  geom_point(data=d_sum,aes(x=group, y=tau),size=2.2) +
  geom_errorbar(data=d_sum,aes(x=group, y=tau,ymin=tau-2*tau_se, ymax=tau+2*tau_se), width=0, position=position_dodge(0.0))+
  geom_line(data=d_sum,aes(x=group, y=tau,group = 1),linewidth=1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black",linewidth=0.5) +
  xlab("Number of shared tissues") +
  ylab(expression("LDSC" ~ tau ~ ", non-coding variants")) +
  scale_x_discrete(labels = c("19", "16-18", "9-15","3-8","1-2")) +
  scale_y_continuous(trans = squish_trans(), breaks = breaks) +
  labs(color = "Trait-Tissue pair") +
  scale_color_manual(name = "Trait-Tissue pair", values = custom_palette) +
  theme_cowplot(9) + 
  guides(color = "none")


outfile="paper_figures/fig4c.pdf"
save_plot(outfile, p0,base_width = 2.3, base_height = 2.1)


d_ldsc=fread("data/ATAC_specificity_LDSC_Mineto.txt")

d_ldsc <- separate(d_ldsc, annot, into = c("tissue", "x", "group"), sep = "_")
d_ldsc <- separate(d_ldsc, x, into = c("covar", "y"), sep = "-") %>% select(-y)

d_data=d_ldsc[d_ldsc$covar=="intensity",]

tissues=sort(unique(d_data$tissue))
traits=sort(unique(d_data$trait))
trait_names=c("MCV","Lymphocyte count","Eosinophil count","Reticulocyte %", "ALT", "Creatinine", "Glucose", "Heel BMD", "Height")
d_map=data.frame(trait=traits,trait_name=trait_names)

d_data=left_join(d_data,d_map)
d_data$trait_cell=paste0(d_data$trait_name,"-",d_data$tissue)

###meta-analysis

d_data$weight <- 1 / (d_data$tau_se)^2

iter=0
for (group_temp in c("GROUP1", "GROUP2", "GROUP3", "GROUP4")){
  iter=iter+1
  d_temp=d_data[d_data$group==group_temp,]
  weighted_mean <- sum(d_temp$weight * d_temp$tau) / sum(d_temp$weight)
  weighted_se <- sqrt(1 / sum(d_temp$weight))
  d_sum_temp=data.frame(group= group_temp,tau=weighted_mean,tau_se=weighted_se,trait_cell="Average")
  if (iter==1){d_sum=d_sum_temp}
  if (iter>1){d_sum=rbind(d_sum,d_sum_temp)}
}

########

d_data$group=factor(d_data$group, levels = rev(c("GROUP1", "GROUP2", "GROUP3", "GROUP4")))
d_sum$group=factor(d_sum$group, levels = rev(c("GROUP1", "GROUP2", "GROUP3", "GROUP4")))

d_all=rbind(d_sum,d_data %>% select(colnames(d_sum)))
d_all$trait_cell <- factor(d_all$trait_cell,levels = c(sort(unique(d_data$trait_cell)),"Average")  )

########

b_point = 0.65
squish_power = 0.1
squish_const = 0.2
squish_trans <- function() {
  trans_new(
    name = "squish",
    transform = function(x) sign(x)*((abs(x)+squish_const)^(squish_power)*squish_const^(squish_power-1)/squish_power - squish_const^(2*squish_power-1)/squish_power),
    inverse = function(x) sign(x)*(((squish_power*abs(x)+squish_const^(2*squish_power-1))/squish_const^(squish_power-1))^(1/squish_power) - squish_const)
  )
}

num_levels <- length(levels(d_all$trait_cell))
custom_palette <- brewer.pal(num_levels - 1, "Set3")
custom_palette <- c(custom_palette, "black")


breaks <- c(-4, -3, -2, -1, 0, 1, 2)
labels <- c(-4, -3, -2, -1, 0, 1, 2)

###

p0 <- ggplot() +
  geom_line(data=d_all,aes(x=group, y=tau,group=trait_cell,color=trait_cell),linewidth=0.5) +
  geom_point(data=d_sum,aes(x=group, y=tau),size=2.2) +
  geom_errorbar(data=d_sum,aes(x=group, y=tau,ymin=tau-2*tau_se, ymax=tau+2*tau_se), width=0, position=position_dodge(0.0))+
  geom_line(data=d_sum,aes(x=group, y=tau,group = 1),linewidth=1) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black",linewidth=0.5) +
  xlab("ATAC-seq peak intensity bin") +
  ylab(expression("LDSC" ~ tau ~ ", non-coding variants")) +
  scale_x_discrete(labels = c("Bin 2", "Bin 3", "Bin 4","Bin 5")) +
  scale_y_continuous(trans = squish_trans(), breaks = breaks) +
  labs(color = "Trait-Tissue pair") +
  scale_color_manual(name = "Trait-Tissue pair", values = custom_palette) +
  theme_cowplot(9) + 
  guides(color = "none")


outfile="paper_figures/figS11.pdf"
save_plot(outfile, p0,base_width = 3.5, base_height = 3)











#########################
### Figure S12 ###
#########################

d_expected=fread("data/zeng_et_al_v3.tsv")
d_cds=read.csv('data/cds_lengths.csv.gz')

d_cds['ensg'] = d_cds$Gene.stable.ID
d_cds = d_cds[c('ensg', 'CDS.Length')]
d_cds <- d_cds[!duplicated(d_cds$ensg),] 

d_joint <- inner_join(d_expected, d_cds, by='ensg')


p1 <- ggplot(d_joint, aes(x=CDS.Length, y=exp_lof)) +
  stat_poly_eq() +
  geom_point(size=2, alpha=0.05) +
  xlab("CDS length in base pairs") +
  theme_cowplot(9) +
  ylab("Expected LoFs") 

outfile=paste0("paper_figures/figS12.pdf")
save_plot(outfile, p1,base_width = 3.5, base_height = 3)




















#########################
### Figure 5A-C ###
#########################



d_indep_traits=fread("data/indep_traits.txt")

d_gamma_unbiased=fread("data/unbiased_gamma_sq_by_trait_MAF1_hgncIDs.txt")
d_gamma_hats=fread("data/gamma_hats_by_trait_MAF1_hgncIDs.txt")
d_z_scores=fread("data/zscores_by_trait_MAF1_hgncIDs.txt")
d_SEs=fread("data/SE_hats_by_trait_MAF1_hgncIDs.txt")
d_s_het=fread("data/zeng_et_al_v3.tsv")

#####

set.seed(0)

traits=colnames(d_z_scores %>% dplyr::select(-hgnc_id))
traits=traits[traits %in% d_indep_traits$trait]

E_bin_size = 100
colnames(d_s_het)[2]="hgnc_id"
d_s_het$exp_lof=replace(d_s_het$exp_lof, is.na(d_s_het$exp_lof), 0) #replace NAs with 0
d_s_het=d_s_het[d_s_het$hgnc_id %in% d_gamma_hats$hgnc_id,]
d_s_het$E_cat=ntile(d_s_het$exp_lof,E_bin_size)

d_gamma2_matrix= d_gamma_unbiased %>% dplyr::select(-hgnc_id)
d_SE_matrix= d_SEs %>% dplyr::select(-hgnc_id) 
d_z_matrix= d_z_scores %>% dplyr::select(-hgnc_id)
d_z_matrix=replace(d_z_matrix, is.na(d_z_matrix), 1) #replace NAs with 1

#compute data points per parameter, i.e. gamma, z, etc
d_test_gamma <- ifelse(abs(d_gamma2_matrix) < 4, 1, 1)
d_test_se <- ifelse(abs(d_SE_matrix) < 4, 1, 1)
d_test_z <- ifelse(abs(d_z_matrix) < 4, 1, 1)

d_data=(d_z_scores %>% dplyr::select(hgnc_id))
d_data$mean_z2=rowMeans(d_z_matrix^2,na.rm = T)
d_data$mean_gamma2=rowMeans(d_gamma2_matrix,na.rm = T)
d_data$mean_SE2=rowMeans(d_SE_matrix^2,na.rm = T)
d_data$count_traits_gamma=rowSums(d_test_gamma,na.rm = T)
d_data$count_traits_se=rowSums(d_test_se,na.rm = T)
d_data$count_traits_z=rowSums(d_test_z,na.rm = T)

d_data=left_join(d_data, d_s_het %>% dplyr::select(hgnc_id,exp_lof,E_cat))
d_data=d_data[!is.na(d_data$E_cat),]

data_sum <- d_data %>%
  group_by(E_cat) %>%
  summarize(mean_gamma_sq = mean(mean_gamma2, na.rm = TRUE),
            mean_se_sq = mean(mean_SE2, na.rm = TRUE),
            mean_z_sq=mean(mean_z2, na.rm = TRUE),
            mean_E = mean(exp_lof, na.rm = TRUE),
            sd_gamma_sq = sd(mean_gamma2, na.rm = TRUE),
            sd_se_sq = sd(mean_SE2, na.rm = TRUE),
            sd_z_sq = sd(mean_z2, na.rm = TRUE),
            count_gamma_sq = sum(count_traits_gamma),
            count_se_sq = sum(count_traits_se),
            count_z_sq = sum(count_traits_z)
  )

data_sum$se_gamma_sq=data_sum$sd_gamma_sq/sqrt(data_sum$count_gamma_sq)
data_sum$se_se_sq=data_sum$sd_se_sq/sqrt(data_sum$count_se_sq)
data_sum$se_z_sq=data_sum$sd_z_sq/sqrt(data_sum$count_z_sq)

#########
plot_scale = 1

span_temp=0.8
color_temp="#006400"
p0 <- ggplot(data_sum, aes(x=mean_E, y=mean_z_sq)) +
  geom_point(size=2,alpha=0.2,color=color_temp) + 
  geom_smooth(method = "loess", span=span_temp,se = F, color=alpha(color_temp, 0.6)) +
  xlab("Mean expected LoFs in bin") +
  ylab(expression("Mean(" * z^2 * ") across traits")) +
  theme_cowplot(9)

outfile="paper_figures/fig5c.pdf"
save_plot(outfile, p0,base_width = 2.25*plot_scale, base_height = 2*plot_scale)

#
color_temp="#4169E1"
p1 <- ggplot(data_sum, aes(x=mean_E, y=mean_gamma_sq)) +
  geom_point(size=2,alpha=0.2,color=color_temp) + 
  geom_smooth(method = "loess", span=span_temp,se = F, color=alpha(color_temp, 0.6)) +
  xlab("Mean expected LoFs in bin") +
  ylab(expression(paste("Mean(" * gamma^2 * ")" ~ "across traits"))) +
  theme_cowplot(9)

outfile="paper_figures/fig5a.pdf"
save_plot(outfile, p1,base_width = 2.25*plot_scale, base_height = 2*plot_scale)

#

color_temp="#DC143C"
p2 <- ggplot(data_sum, aes(x=mean_E, y=mean_se_sq)) +
  geom_point(size=2,alpha=0.2,color=color_temp) + 
  geom_smooth(method = "loess", span=span_temp,se = F, color=alpha(color_temp, 0.6)) +
  xlab("Mean expected LoFs in bin") +
  ylab(expression(paste("Mean SE(" * hat(gamma) * ")"^2 ~ "across traits"))) +
  theme_cowplot(9)

outfile="paper_figures/fig5b.pdf"
save_plot(outfile, p2,base_width = 2.25*plot_scale, base_height = 2*plot_scale)













#########################
### Figure S13 ###
#########################

trajectories <- suppressWarnings(melt(as.data.table(t(read.table('data/wright_fisher_trajectories.txt')/20000))))

trajectories$time <- rep(1:1000, 10000)


p1 <- ggplot(trajectories, aes(x=time, y=value, group=variable)) + 
  rasterize(geom_line(color=alpha('black', 0.2)), dpi=600) + 
  xlab("Generations since mutation arose") +
  ylab("Frequency") +
  scale_x_continuous(breaks=c(0, 250, 500, 750, 1000),labels=c(0, 250, 500, 750, 'Present Day'), expand=expansion(add=c(50, 100)))

outfile="paper_figures/figS13.pdf"
save_plot(outfile, p1, base_width = 4, base_height = 2)



















#########################
### Figure 5D ###
#########################


sims <- data.frame(read.table('data/simulated_realized_heritability.txt'))

sims <- sims[sims$V2 > 0, ]

sims$V1 <- sims$V1 / max(sims$V1)
sims$V2 <- sims$V2 / max(sims$V2)
sims$MAF <- sims$V2 / sims$V1
sims$MAF<- sims$MAF* 0.482355 / max(sims$MAF)
sims$psi <- rep(1, dim(sims)[1])

p1 <- ggplot(sims, aes(x=V1, y=V2, colour=MAF)) + 
  geom_point() +
  scale_color_viridis(direction = -1) +
  xlab("Scaled squared effect on study trait") +
  ylab("Realized heritability (relative)")

outfile="paper_figures/fig5d.pdf"
save_plot(outfile, p1, base_width = 3.5*1.0, base_height = 2.25*1.2)











#########################
### Figure 5F-H ###
#########################


d_indep_traits=fread("data/indep.phenos_irnt.txt",header = F)
colnames(d_indep_traits)[1]="pheno"

#GWAS summary stats for all Backman traits
d_pval=fread("data/pval_by_trait.txt")
shared_traits=intersect(d_indep_traits$pheno,colnames(d_pval))
plei <- rowSums((d_pval %>% select(shared_traits) < 5e-8)[, -1])
d_pval$plei <- plei
d_alpha=fread("data/alpha_hats_by_trait.txt.gz")
d_SE=fread("data/alpha_SEs_by_trait.txt.gz")

d_maf <- fread('data/SNP.maf.txt.gz')



d_pruned_hits=fread("data/gwas_hits.blocks.txt")

d_pruned_hits <- merge(d_maf, d_pruned_hits, by.x = 'SNP', by.y = 'SNP', all.x = F, all.y = T)

d_hit_counts=aggregate(d_pruned_hits$SNP,list(d_pruned_hits$pheno),length)
colnames(d_hit_counts)=c("pheno","hit_count")
d_hit_counts=d_hit_counts %>% arrange(-hit_count)
d_hit_counts=d_hit_counts[d_hit_counts$pheno %in% shared_traits,] 
rownames(d_hit_counts)=1:nrow(d_hit_counts)

#####
#compute pleiotropy per hit-trait pair (for how many other traits the SNP is a hit)
d_stats=data.frame(hit_cat=as.numeric(),one=as.numeric(),two=as.numeric(),three=as.numeric(),four=as.numeric(), total=as.numeric(), trait =as.character())

#18 traits with >100 hits
for (i in 1:18){
  print(i)
  pheno_temp=d_hit_counts$pheno[i]
  d_hits_temp=d_pruned_hits[d_pruned_hits$pheno==pheno_temp,]
  d_hits_temp$pcat=ntile(d_hits_temp$pval,4) #4 rank bins
  
  #loop over rank bins
  for (k in 1:4){
    hits=d_hits_temp[d_hits_temp$pcat==k,]
    d_pval_temp=d_pval[d_pval$SNP %in% hits$SNP,]
    d_pval_indep_temp=d_pval_temp %>% select(shared_traits)
    
    dx=as.matrix(d_pval_indep_temp)
    gwas_thresh= 5e-8
    dx[dx <= gwas_thresh] <- 2
    dx[dx != 2] <- 0
    dx[dx == 2] <- 1
    
    d_plei_temp=data.frame(SNP=d_pval_temp$SNP,plei=rowSums(dx))
    d_plei_temp=d_plei_temp[d_plei_temp$plei>0,] %>% arrange(-plei)
    x5=mean(d_plei_temp$plei)
    x6=mean(hits$maf_neale, na.rm=TRUE)
    x=c(k, x5, x6)
    names(x)=c("hit_cat", "mean_traits", "mean_maf")
    
    dx=data.frame(t(x))
    dx$trait=pheno_temp
    d_stats=rbind(d_stats,dx)
  }
  
}

#####
#average across traits

d_stats %>% 
  group_by(hit_cat) %>% 
  summarize(mean_traits=mean(mean_traits), mean_maf=mean(mean_maf)) -> d_sum

theory <- read.csv('data/p_val_ranking_sims.csv')[,c('rep', 'p_bin_1_mean_num_hits', 'p_bin_2_mean_num_hits', 'p_bin_3_mean_num_hits', 'p_bin_4_mean_num_hits')]
theory <- melt(as.data.table(theory), id.vars='rep')
theory$p_bin = rep(NA, dim(theory)[1])
theory$p_bin[grepl('p_bin_1', theory$variable)] = 1
theory$p_bin[grepl('p_bin_2', theory$variable)] = 2
theory$p_bin[grepl('p_bin_3', theory$variable)] = 3
theory$p_bin[grepl('p_bin_4', theory$variable)] = 4



theory_spec <- read.csv('data/p_val_ranking_sims.csv')[,c('rep', 'p_bin_1_specificity', 'p_bin_2_specificity', 'p_bin_3_specificity', 'p_bin_4_specificity')]
theory_spec <- melt(as.data.table(theory_spec), id.vars='rep')
theory_spec$p_bin = rep(NA, dim(theory_spec)[1])
theory_spec$p_bin[grepl('p_bin_1', theory_spec$variable)] = 1
theory_spec$p_bin[grepl('p_bin_2', theory_spec$variable)] = 2
theory_spec$p_bin[grepl('p_bin_3', theory_spec$variable)] = 3
theory_spec$p_bin[grepl('p_bin_4', theory_spec$variable)] = 4


theory_freq <- read.csv('data/p_val_ranking_sims.csv')[,c('rep', 'p_bin_1_freq', 'p_bin_2_freq', 'p_bin_3_freq', 'p_bin_4_freq')]
theory_freq <- melt(as.data.table(theory_freq), id.vars='rep')
theory_freq$p_bin = rep(NA, dim(theory_freq)[1])
theory_freq$p_bin[grepl('p_bin_1', theory_freq$variable)] = 1
theory_freq$p_bin[grepl('p_bin_2', theory_freq$variable)] = 2
theory_freq$p_bin[grepl('p_bin_3', theory_freq$variable)] = 3
theory_freq$p_bin[grepl('p_bin_4', theory_freq$variable)] = 4


#####

p1<-ggplot(d_sum, aes(x = hit_cat, y = mean_traits)) +
  geom_line(data=theory, aes(x=p_bin, y=value, group=rep), color=alpha("orange", 0.1)) + 
  geom_point(colour = "black") +
  xlab("p-value rank") +
  ylab("Mean number of\nsignificant traits per hit")

outfile=paste0("paper_figures/fig5h.pdf")
save_plot(outfile, p1,base_width = 2, base_height = 2)


d_mean_maf <- d_stats %>% group_by(hit_cat) %>% summarise(mean_maf=mean(mean_maf) / mean(d_stats$mean_maf))

p2 <- ggplot(d_mean_maf, aes(x = hit_cat, y = mean_maf)) +
  geom_line(data=theory_freq %>% group_by(rep) %>% mutate(value=value/mean(value)), aes(x=p_bin, y=value, group=rep), color=alpha("orange", 0.1)) +
  geom_point() +
  xlab("p-value rank") +
  ylab("Mean MAF relative\nto overall mean MAF")

outfile=paste0("paper_figures/fig5f.pdf")
save_plot(outfile, p2,base_width = 2, base_height = 2)


p3 <- ggplot(theory_spec, aes(x = p_bin, y = value, group=rep)) +
  geom_line(color=alpha("orange", 0.1)) +
  xlab("p-value rank") +
  ylab("Mean specificity")

outfile=paste0("paper_figures/fig5g.pdf")
save_plot(outfile, p3,base_width = 1.85, base_height = 2)









#########################
### Figure S14-17 ###
#########################


make_plots <- function(N, p, sig_p, scale_factor){
  
  fname <- paste0('data/p_val_ranking_sims_N_', N, '_p_', p, '_sig_p_', sig_p, '_scale_factor_', scale_factor, '.csv')
  
  theory <- read.csv(fname)[,c('rep', 'p_bin_1_mean_num_hits', 'p_bin_2_mean_num_hits', 'p_bin_3_mean_num_hits', 'p_bin_4_mean_num_hits')]
  theory <- melt(as.data.table(theory), id.vars='rep')
  theory$p_bin = rep(NA, dim(theory)[1])
  theory$p_bin[grepl('p_bin_1', theory$variable)] = 1
  theory$p_bin[grepl('p_bin_2', theory$variable)] = 2
  theory$p_bin[grepl('p_bin_3', theory$variable)] = 3
  theory$p_bin[grepl('p_bin_4', theory$variable)] = 4
  
  
  
  theory_spec <- read.csv(fname)[,c('rep', 'p_bin_1_specificity', 'p_bin_2_specificity', 'p_bin_3_specificity', 'p_bin_4_specificity')]
  theory_spec <- melt(as.data.table(theory_spec), id.vars='rep')
  theory_spec$p_bin = rep(NA, dim(theory_spec)[1])
  theory_spec$p_bin[grepl('p_bin_1', theory_spec$variable)] = 1
  theory_spec$p_bin[grepl('p_bin_2', theory_spec$variable)] = 2
  theory_spec$p_bin[grepl('p_bin_3', theory_spec$variable)] = 3
  theory_spec$p_bin[grepl('p_bin_4', theory_spec$variable)] = 4
  
  
  theory_freq <- read.csv(fname)[,c('rep', 'p_bin_1_freq', 'p_bin_2_freq', 'p_bin_3_freq', 'p_bin_4_freq')]
  theory_freq <- melt(as.data.table(theory_freq), id.vars='rep')
  theory_freq$p_bin = rep(NA, dim(theory_freq)[1])
  theory_freq$p_bin[grepl('p_bin_1', theory_freq$variable)] = 1
  theory_freq$p_bin[grepl('p_bin_2', theory_freq$variable)] = 2
  theory_freq$p_bin[grepl('p_bin_3', theory_freq$variable)] = 3
  theory_freq$p_bin[grepl('p_bin_4', theory_freq$variable)] = 4
  
  
  #####
  
  p1<-ggplot(theory, aes(x = p_bin, y = value, group=rep)) +
    geom_line(color=alpha("orange", 0.1)) + 
    geom_line(aes(x = p_bin, y = value, group=NA), data=theory%>%group_by(p_bin)%>%summarize(value=mean(value))) +
    xlab("p-value rank") +
    ylab("Mean number of\nsignificant traits per hit")
  
  p2 <- ggplot(theory_freq %>% group_by(rep) %>% mutate(value=value/mean(value)), aes(x=p_bin, y=value, group=rep)) +
    geom_line(color=alpha("orange", 0.1)) +
    geom_line(aes(x = p_bin, y = value, group=NA), data=theory_freq%>% group_by(rep) %>% mutate(value=value/mean(value)) %>% ungroup() %>% group_by(p_bin)%>%summarize(value=mean(value))) +
    xlab("p-value rank") +
    ylab("Mean MAF relative\nto overall mean MAF")
  
  p3 <- ggplot(theory_spec, aes(x = p_bin, y = value, group=rep)) +
    geom_line(color=alpha("orange", 0.1)) +
    geom_line(aes(x = p_bin, y = value, group=NA), data=theory_spec%>%group_by(p_bin)%>%summarize(value=mean(value))) +
    xlab("p-value rank") +
    ylab("Mean specificity")
  
  fname <- paste0('data/p_val_ranking_sims_hist_N_', N, '_p_', p, '_sig_p_', sig_p, '_scale_factor_', scale_factor, '.csv')
  theory <- read.csv(fname)
  p4 <- ggplot(theory, aes(x=specificity)) +
    geom_histogram() +
    xlab("Specificity") + 
    ylab("Num. seg. sites")
  
  to_return = list()
  to_return[[1]] = p1
  to_return[[2]] = p2
  to_return[[3]] = p3
  to_return[[4]] = p4
  return(to_return)
}


freq_plots <- list()
spec_plots <- list()
plei_plots <- list()
hist_plots <- list()
for(i in 1:4){
  N <- c('5000000.0', '10000000', '20000000', '100000000')[i]
  these_plots <- make_plots(N, '0.5', '1e-05', '0.33')
  freq_plots[[i]] <- these_plots[[2]]
  spec_plots[[i]] <- these_plots[[3]]
  plei_plots[[i]] <- these_plots[[1]]
  hist_plots[[i]] <- these_plots[[4]]
}
p_full <- plot_grid(freq_plots[[1]], freq_plots[[2]], freq_plots[[3]], freq_plots[[4]],
                    spec_plots[[1]], spec_plots[[2]], spec_plots[[3]], spec_plots[[4]],
                    plei_plots[[1]], plei_plots[[2]], plei_plots[[3]], plei_plots[[4]],
                    hist_plots[[1]], hist_plots[[2]], hist_plots[[3]], hist_plots[[4]],
                    nrow=4, ncol=4, byrow=FALSE)

save_plot('paper_figures/figS14.pdf', p_full,base_width = 2.25*4, base_height = 2*4)




freq_plots <- list()
spec_plots <- list()
plei_plots <- list()
hist_plots <- list()
for(i in 1:5){
  p <- c('0.1', '0.25', '0.5', '0.75', '1.05')[i]
  these_plots <- make_plots('10000000', p, '1e-05', '0.33')
  freq_plots[[i]] <- these_plots[[2]]
  spec_plots[[i]] <- these_plots[[3]]
  plei_plots[[i]] <- these_plots[[1]]
  hist_plots[[i]] <- these_plots[[4]]
}
p_full <- plot_grid(freq_plots[[1]], freq_plots[[2]], freq_plots[[3]], freq_plots[[4]], freq_plots[[5]],
                    spec_plots[[1]], spec_plots[[2]], spec_plots[[3]], spec_plots[[4]], spec_plots[[5]],
                    plei_plots[[1]], plei_plots[[2]], plei_plots[[3]], plei_plots[[4]], plei_plots[[5]],
                    hist_plots[[1]], hist_plots[[2]], hist_plots[[3]], hist_plots[[4]], hist_plots[[5]],
                    nrow=5, ncol=4, byrow=FALSE)

save_plot('paper_figures/figS16.pdf', p_full,base_width = 2.25*4, base_height = 2*5)



freq_plots <- list()
spec_plots <- list()
plei_plots <- list()
hist_plots <- list()
for(i in 1:5){
  sig_p <- c('0.001', '0.0001', '1e-05', '1e-06', '1e-07')[i]
  these_plots <- make_plots('10000000', '0.5', sig_p, '0.33')
  freq_plots[[i]] <- these_plots[[2]]
  spec_plots[[i]] <- these_plots[[3]]
  plei_plots[[i]] <- these_plots[[1]]
  hist_plots[[i]] <- these_plots[[4]]
}
p_full <- plot_grid(freq_plots[[1]], freq_plots[[2]], freq_plots[[3]], freq_plots[[4]], freq_plots[[5]],
                    spec_plots[[1]], spec_plots[[2]], spec_plots[[3]], spec_plots[[4]], spec_plots[[5]],
                    plei_plots[[1]], plei_plots[[2]], plei_plots[[3]], plei_plots[[4]], plei_plots[[5]],
                    hist_plots[[1]], hist_plots[[2]], hist_plots[[3]], hist_plots[[4]], hist_plots[[5]],
                    nrow=5, ncol=4, byrow=FALSE)

save_plot('paper_figures/figS15.pdf', p_full,base_width = 2.25*4, base_height = 2*5)



freq_plots <- list()
spec_plots <- list()
plei_plots <- list()
hist_plots <- list()
for(i in 1:5){
  scale_factor <- c('0.0033', '0.033', '0.33', '3.3000000000000003', '33.0')[i]
  these_plots <- make_plots('10000000', '0.5', '1e-05', scale_factor)
  freq_plots[[i]] <- these_plots[[2]]
  spec_plots[[i]] <- these_plots[[3]]
  plei_plots[[i]] <- these_plots[[1]]
  hist_plots[[i]] <- these_plots[[4]]
}
p_full <- plot_grid(freq_plots[[1]], freq_plots[[2]], freq_plots[[3]], freq_plots[[4]], freq_plots[[5]],
                    spec_plots[[1]], spec_plots[[2]], spec_plots[[3]], spec_plots[[4]], spec_plots[[5]],
                    plei_plots[[1]], plei_plots[[2]], plei_plots[[3]], plei_plots[[4]], plei_plots[[5]],
                    hist_plots[[1]], hist_plots[[2]], hist_plots[[3]], hist_plots[[4]], hist_plots[[5]],
                    nrow=5, ncol=4, byrow=FALSE)

save_plot('paper_figures/figS17.pdf', p_full,base_width = 2.25*4, base_height = 2*5)










#########################
### Figure 6A ###
#########################


theta <- 0.0001
expected_freq_var <- function(ns){
  to_return <- exp(
    BAS::hypergeometric1F1(1/2, theta + 3/2, ns, log=TRUE) +
      - log(8) +
      - log(theta + 1/2) +
      - BAS::hypergeometric1F1(1/2, theta+1/2, ns, log=TRUE)
  )
  return(to_return)
}

# xs <- seq(0, 5, length.out=500)
xs <-exp(seq(log(0.05), 10, length.out=500))

to_plot <- matrix(nrow=20, ncol=500)

psis <- seq(0.05, 1, length.out=20)
for(k in 1:20){
  these_xs <- xs/psis[k]
  these_xs[these_xs > 80] <- 80
  these_mods = these_xs * psis[k]
  to_plot[k, ] <- 8 * these_mods * sapply(these_xs, expected_freq_var)
}

to_plot <- suppressWarnings(melt(as.data.table(t(to_plot))))
to_plot$time <- rep(xs/5, 20)
to_plot$psi <- rep(psis, each=500)

p1 <- ggplot(to_plot, aes(x=time, y=value, group=psi)) + 
  geom_line(aes(color=psi)) +
  scale_color_viridis(direction = -1) +
  # xlab("Scaled squared effect on study trait") +
  xlab("") +
  scale_x_continuous(trans='log') +
  # scale_y_continuous(trans='log') +
  ylab(expression("Expected contribution to " ~ h^2)) +
  theme(axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank()) +
  labs(colour = expression(Psi["v"]))

outfile="paper_figures/fig6a.pdf"
save_plot(outfile, p1, base_width = 3.0, base_height = 1.945)












#########################
### Figure 6D ###
#########################



df <- as.data.table(read.csv('data/amm_results.csv'))

df_means <- df %>% group_by(s_het) %>% summarise(ivw_h2=sum(h2_mean / h2_se^2 / (sum(1/h2_se^2))), raw_mean=mean(h2_mean), se=sqrt(1/(sum(1/h2_se^2))))

df_means$ivw_h2 <- df_means$ivw_h2 / mean(df_means$ivw_h2)

p1 <- ggplot(df_means, aes(x=s_het, y=ivw_h2)) +
  geom_point(size=2,alpha=0.2, color='purple') +
  geom_smooth(method = "loess", span=span_temp,se = F, color = alpha('purple', 0.6)) +
  ylab('Heritability enrichment') +
  xlab(substitute(paste("Mean s"[het], " in bin"))) +
  scale_x_continuous(trans='log', breaks=c(1e-4, 1e-3, 1e-2, 1e-1, 1))


outfile="paper_figures/fig6d.pdf"
save_plot(outfile, p1,base_width = 2.5, base_height = 2)










#########################
### Figure S18 ###
#########################




d_data=fread("data/logistic_reg_gwas_hits_by_s_het_remove_pa.txt")

genefile="data/zeng_et_al_v3.tsv"
d_gene=fread(genefile)
d_gene=d_gene %>% dplyr::select(hgnc,post_mean)
colnames(d_gene)[1]="gene"
d_gene$s_het_cat=ntile(d_gene$post_mean,100)

data_sum <- d_gene %>%
  group_by(s_het_cat) %>%
  summarize(mean_s_het = mean(post_mean, na.rm = TRUE))

d_data$mean_s_het=data_sum$mean_s_het

#####

span_temp=0.8
line_alpha=0.8
point_alpha=0.2

color_temp="purple"
p1 <- ggplot(d_data, aes(x=mean_s_het, y=coeff)) +
  geom_point(size=2,alpha=point_alpha,color=color_temp) + 
  geom_smooth(method = "loess", span=span_temp,se = F, color = alpha(color_temp, line_alpha)) +
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = trans_format("log10", math_format(10^.x))) +
  xlab(substitute(paste("Mean s"[het], " in bin"))) +
  ylab("Log odds of GWAS hit") +
  theme_cowplot(9)

outfile=("paper_figures/figS18.pdf")
save_plot(outfile, p1,base_width = 3.5, base_height = 3)






#########################
### Figure S19 ###
#########################


df <- read_csv('data/gwas_num_hits_correlation_with_gamma_sq.csv')

span_temp=0.8
line_alpha=0.8
point_alpha=0.2


p0 <- ggplot(df, aes(x=num_hits, y=all_corrs)) +
  geom_hline(yintercept=0, color='red', linetype=2) +
  geom_point(size=2, alpha=point_alpha) +
  xlab('Total number of independent GWAS hits') +
  ylab(expression("Correlaton of number of hits with "~widehat(gamma^2))) +
  ggtitle("All genes")



p1 <- ggplot(df, aes(x=num_hits, y=corr_just_hits)) +
  geom_hline(yintercept=0, color='red', linetype=2) +
  geom_point(size=2, alpha=point_alpha) +
  xlab('Total number of independent GWAS hits') +
  ylab("Correlaton of number of hits with "~widehat(gamma^2)) +
  ggtitle("Genes with at least one hit")


outfile=("paper_figures/figS19.pdf")
pfull <- plot_grid(p0, p1, ncol=2, nrow=1)

save_plot(outfile, pfull, base_width = 9, base_height = 4)
