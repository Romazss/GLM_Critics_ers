library(GLMsData) # Load the GLMsData package
data(lungcap)
head(lungcap)

str(lungcap)

lungcap$Smoke <- factor(lungcap$Smoke,
                        levels=c(0, 1), # The values of Smoke
                        labels=c("Non-smoker","Smoker")) # The labels

summary(lungcap$Smoke)


plot( FEV ~ Ht, data=lungcap, main="FEV vs height",
      xlab="Height (in inches)", ylab="FEV (in L)",
      las=1, ylim=c(0, 6) )
plot( FEV ~ Gender, data=lungcap,
        main="FEV vs gender", ylab="FEV (in L)",
        las=1, ylim=c(0, 6))
plot( FEV ~ Smoke, data=lungcap, main="FEV vs Smoking status",
        ylab="FEV (in L)", xlab="Smoking status",
        las=1, ylim=c(0, 6))


plot( FEV ~ Age,
      data=subset(lungcap, Smoke=="Smoker"), # Only select smokers
      main="FEV vs age\nfor smokers", # \n means `new line'
      ylab="FEV (in L)", xlab="Age (in years)",
      ylim=c(0, 6), xlim=c(0, 20), las=1)

plot( FEV ~ Age,
      data=subset(lungcap, Smoke=="Non-smoker"), # Only select non-smokers
      main="FEV vs age\nfor non-smokers",
      ylab="FEV (in L)", xlab="Age (in years)",
      ylim=c(0, 6), xlim=c(0, 20), las=1)

plot( FEV ~ Ht, data=subset(lungcap, Smoke=="Smoker"),
        main="FEV vs height\nfor smokers",
        ylab="FEV (in L)", xlab="Height (in inches)",
        xlim=c(45, 75), ylim=c(0, 6), las=1)

plot( FEV ~ Ht, data=subset(lungcap, Smoke=="Non-smoker"),
        main="FEV vs height\nfor non-smokers",
        ylab="FEV (in L)", xlab="Height (in inches)",
        xlim=c(45, 75), ylim=c(0, 6), las=1)

################# PARA SYLLABUS #############
# Ruta donde quieres guardar los PDFs
ruta <- "C:/Users/Ernesto San Martin/Dropbox/MisDocumentosNoteBook/Docencia/GLM_Criticis/GLM2026/"

# =========================
# PDF 1: relaciones generales
# =========================
pdf(paste0(ruta, "lungcap_general.pdf"), width=8, height=10)

par(mfrow=c(3,1))

plot(FEV ~ Ht, data=lungcap,
     main="FEV vs estatura",
     xlab="Estatura (pulgadas)",
     ylab="FEV (litros)",
     las=1, ylim=c(0,6))

plot(FEV ~ Gender, data=lungcap,
     main="FEV vs género",
     xlab="Género",
     ylab="FEV (litros)",
     las=1, ylim=c(0,6))

plot(FEV ~ Smoke, data=lungcap,
     main="FEV vs hábito de fumar",
     xlab="Estado de fumador",
     ylab="FEV (litros)",
     las=1, ylim=c(0,6))

dev.off()


# =========================
# PDF 2: edad por hábito de fumar
# =========================
pdf(paste0(ruta, "lungcap_edad.pdf"), width=8, height=8)

par(mfrow=c(2,1))

plot(FEV ~ Age,
     data=subset(lungcap, Smoke=="Smoker"),
     main="FEV vs edad (fumadores)",
     xlab="Edad (años)",
     ylab="FEV (litros)",
     ylim=c(0,6), xlim=c(0,20), las=1)

plot(FEV ~ Age,
     data=subset(lungcap, Smoke=="Non-smoker"),
     main="FEV vs edad (no fumadores)",
     xlab="Edad (años)",
     ylab="FEV (litros)",
     ylim=c(0,6), xlim=c(0,20), las=1)

dev.off()


# =========================
# PDF 3: estatura por hábito de fumar
# =========================
pdf(paste0(ruta, "lungcap_estatura.pdf"), width=8, height=8)

par(mfrow=c(2,1))

plot(FEV ~ Ht,
     data=subset(lungcap, Smoke=="Smoker"),
     main="FEV vs estatura (fumadores)",
     xlab="Estatura (pulgadas)",
     ylab="FEV (litros)",
     xlim=c(45,75), ylim=c(0,6), las=1)

plot(FEV ~ Ht,
     data=subset(lungcap, Smoke=="Non-smoker"),
     main="FEV vs estatura (no fumadores)",
     xlab="Estatura (pulgadas)",
     ylab="FEV (litros)",
     xlim=c(45,75), ylim=c(0,6), las=1)

dev.off()
