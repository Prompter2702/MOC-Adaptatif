# =========================
# Configuration utilisateur
# =========================

FC      = gfortran
SRCDIR   = src
BUILDDIR = build

MODSRC = $(SRCDIR)/flgbcd.for $(SRCDIR)/flginc.for $(SRCDIR)/flgoct.for \
		$(SRCDIR)/flgsnq.for $(SRCDIR)/flxcom.for $(SRCDIR)/snqhrm.for \
		$(SRCDIR)/sweep8one.for $(SRCDIR)/srccorr.for
SRC        = $(wildcard  $(SRCDIR)/*.for)
SRCSANSMOD = $(filter-out $(MODSRC), $(SRC))

OBJ = $(patsubst $(SRCDIR)/%.for, $(BUILDDIR)/%.o, $(MODSRC)) \
      $(patsubst $(SRCDIR)/%.for, $(BUILDDIR)/%.o, $(SRCSANSMOD))
		   
EXE     = sweeprecur
MODDIR  = moddir

# =========================
# Compilation options
# =========================

FFLAGS_DEBUG   = -O0 -Wall -fcheck=all -finit-real=snan -g -J$(MODDIR) -I$(MODDIR)
FFLAGS_FAST    = -O2 -march=native -funroll-loops -Wall -J$(MODDIR) -I$(MODDIR)
FFLAGS_CURRENT = $(FFLAGS_DEBUG)
# FFLAGS_CURRENT = $(FFLAGS_FAST)

# =========================
# Règles
# =========================

$(EXE): $(OBJ)
	$(FC) $(FFLAGS_CURRENT) -o $@ $^

$(BUILDDIR)/%.o: $(SRCDIR)/%.for | $(BUILDDIR) $(MODDIR) 
	$(FC) $(FFLAGS_CURRENT) -c $< -o $@

$(BUILDDIR):
	mkdir -p $(BUILDDIR)

# Crée le dossier des modules si absent
$(MODDIR):
	mkdir -p $(MODDIR)

# =========================
# Cibles utiles
# =========================

.PHONY: clean debug release run

clean:
	rm -f $(BUILDDIR)/*.o $(MODDIR)/*.mod $(EXE)
	rm -rf $(MODDIR) $(BUILDDIR)
	rm fort.9

debug:
	$(MAKE) FFLAGS_CURRENT="$(FFLAGS_DEBUG)"

release:
	$(MAKE) FFLAGS_CURRENT="$(FFLAGS_FAST)"

run: $(EXE)
	./$(EXE)


# # Nom de l'exécutable final
# EXEC = sweeprecur

# # Compilateur Fortran
# FC = gfortran

# # Options de compilation
# # FFLAGS = -fcheck=all -O3 -Wall -J$(MODDIR)  # -Jmod place les .mod dans mod/
# FFLAGS = -O2 -march=native -funroll-loops -Wall -J$(MODDIR)  # -Jmod place les .mod dans mod/

# # Dossiers
# SRCDIR = src
# BUILDDIR = build
# MODDIR = mod  # Dossier pour les .mod

# # Liste des fichiers sources
# SRCS = $(wildcard $(SRCDIR)/*.for)

# # Fichiers objets générés après compilation
# OBJS = $(patsubst $(SRCDIR)/%.for, $(BUILDDIR)/%.o, $(SRCS))

# # Règle principale
# $(EXEC): $(OBJS) | $(MODDIR)
# 	$(FC) $(FFLAGS) -o $(EXEC) $(OBJS)

# # Compilation de chaque fichier .for en .o
# $(BUILDDIR)/%.o: $(SRCDIR)/%.for | $(BUILDDIR)
# 	$(FC) $(FFLAGS) -c $< -o $@

# # Création des dossiers si nécessaire
# $(BUILDDIR):
# 	mkdir -p $(BUILDDIR)

# $(MODDIR):
# 	mkdir -p $(MODDIR)

# # Nettoyage
# clean:
# 	rm -rf $(BUILDDIR) $(EXEC) $(MODDIR)
