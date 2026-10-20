import imageio.v2 as imageio
import os

# Chemin vers le dossier contenant les images PNG
chemin_dossier_images = "/Users/emmagurcan/Documents/France/ENS/M1/stageM1/TOEddies/TOEddiesMSM40/CodesArtemis/Figure bathy+ADT/section5"

# Liste pour stocker les chemins des images
chemins_images = []
print("Fichiers trouvés :", os.listdir(chemin_dossier_images))


# Parcours du dossier et ajout des chemins des images à la liste
for fichier in sorted(os.listdir(chemin_dossier_images)):
    if fichier.endswith(".png"):
        chemin_image = os.path.join(chemin_dossier_images, fichier)
        print("Ajout de l'image :", chemin_image)
        chemins_images.append(chemin_image)


# Création de l'animation
with imageio.get_writer('animation.gif', mode='I', duration=0.5) as writer:
    for chemin_image in chemins_images:
        image = imageio.imread(chemin_image)
        writer.append_data(image)

print("L'animation a été créée avec succès!")

# import imageio
# import os

# # Chemin vers le dossier contenant les images PNG
# chemin_dossier_images = "/Users/emmagurcan/Documents/France/ENS/M1/stageM1/from_julian/TOEeddies/TOEddiesMSM40/CodesArtemis/Figure bathy+ADT_eg"

# # Liste pour stocker les chemins des images
# chemins_images = []

# # Parcours du dossier et ajout des chemins des images à la liste
# for fichier in os.listdir(chemin_dossier_images):
#     if fichier.endswith(".png"):
#         chemin_image = os.path.join(chemin_dossier_images, fichier)
#         chemins_images.append(chemin_image)

# # Création de l'animation
# with imageio.get_writer('animation.gif', mode='I', duration=0.5) as writer:
#     for chemin_image in chemins_images:
#         image = imageio.imread(chemin_image)
#         writer.append_data(image)

# print("L'animation a été créée avec succès!")