library(ggplot2)
library(ggpubr)
library(hexSticker)

df<-data.frame(snps=c(53000,39000,32000,28000,26000,22000,20000,18000,12000,8000,2000),
               cutoff=c(.3,.5,.6,.65,.7,.75,.8,.85,.9,.95,1))
plot1<-ggplot2::ggplot(df, ggplot2::aes(x=cutoff)) +
  ggplot2::scale_x_continuous(breaks=c(.3,.5,.6,.65,.7,.75,.8,.85,.9,.95,1))+
  ggplot2::geom_point(ggplot2::aes(y=snps)) +
  ggplot2::theme_bw() + theme_transparent()+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), axis.text.x=element_blank(), axis.text.y=element_blank())+
  ggplot2::labs(x = "SNP completeness cutoff", y = "total SNPs retained")+
  ggplot2::geom_vline(xintercept = .85, color = "red")

plot2<-ggplot2::ggplot(df, ggplot2::aes(x=cutoff)) +
  ggplot2::scale_x_continuous(breaks=c(.3,.5,.6,.65,.7,.75,.8,.85,.9,.95,1))+
  ggplot2::geom_point(ggplot2::aes(y=snps)) +
  ggplot2::theme_bw() + theme_transparent()+
  theme(axis.title.x=element_blank(), axis.title.y=element_blank(), text = element_text(size=10),
        axis.text.x = element_text(angle=45, vjust = 1))+
  #ggplot2::labs(x = "SNP completeness cutoff", y = "total SNPs retained")+
  ggplot2::geom_vline(xintercept = .85, color = "red")


sticker(plot1, package="SNPfiltR", p_size=30, s_x=1, s_y=.85, s_width=1.55, s_height=.88,
        h_fill="cornsilk", h_color="dodgerblue", p_color="dodgerblue", filename="~/Downloads/fggplot2.png")

sticker(plot1, package="SNPfiltR", p_size=5.5, s_x=1, s_y=.85, s_width=1.55, s_height=.88,
        h_fill="#fadadd", h_color="#daedfa", p_color="#dafae7", filename="~/Downloads/fggplot2.png")

sticker(plot1, package="SNPfiltR", p_size=5.5, s_x=1, s_y=.85, s_width=1.55, s_height=.88,
        h_fill="#daedfa", h_color="#dafae7", p_color="#fadadd", filename="~/Downloads/fggplot2.png")

sticker(plot2, package="SNPfiltR", p_size=30, s_x=1, s_y=.85, s_width=1.6, s_height=.9,
        h_fill="#fdf0f2", h_color="#f1979f", p_color="#f7c4c8", filename="~/Downloads/fggplot2.png")

sticker(plot2, package="SNPfiltR", p_size=30, s_x=1, s_y=.85, s_width=1.6, s_height=.9,
        h_fill="#fdf0f2", h_color="#f1979f", p_color="#f7c4c8", filename="~/Desktop/SNPfiltR/man/figures/logo.png")


