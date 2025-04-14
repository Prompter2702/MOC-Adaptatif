# Nom de l'exécutable final
EXEC = sweeprecur

# Compilateur Fortran
FC = gfortran

# Options de compilation
# FFLAGS = -fcheck=all -g -O0 -Wall -J$(MODDIR)  # -Jmod place les .mod dans mod/
FFLAGS = -fcheck=all -g -O0 -Wall -J$(MODDIR)  # -Jmod place les .mod dans mod/

# Dossiers
SRCDIR = src
BUILDDIR = build
MODDIR = mod  # Dossier pour les .mod

# Liste des fichiers sources
SRCS = $(wildcard $(SRCDIR)/*.for)

# Fichiers objets générés après compilation
OBJS = $(patsubst $(SRCDIR)/%.for, $(BUILDDIR)/%.o, $(SRCS))

# Règle principale
$(EXEC): $(OBJS) | $(MODDIR)
	$(FC) $(FFLAGS) -o $(EXEC) $(OBJS)

# Compilation de chaque fichier .for en .o
$(BUILDDIR)/%.o: $(SRCDIR)/%.for | $(BUILDDIR)
	$(FC) $(FFLAGS) -c $< -o $@

# Création des dossiers si nécessaire
$(BUILDDIR):
	mkdir -p $(BUILDDIR)

$(MODDIR):
	mkdir -p $(MODDIR)

# Nettoyage
clean:
	rm -rf $(BUILDDIR) $(EXEC) $(MODDIR)
