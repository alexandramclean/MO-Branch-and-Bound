
/* Fichier effectuant un appel à la fonction combo, avec pour but de rendre plus simple un wrapper pour Julia */

#include "combo.c"

stype comboJulia(item *prob, int nbItem, stype capa, stype maxZ)
{
  return combo(prob,&(prob[nbItem-1]),capa,0,maxZ,1,1);
}

stype comboJuliaMoche(itype *probp, itype *probw, boolean *probx, ntype *probi, int nbItem, stype capa, stype maxZ)
{
    item *prob;
    int i;

    prob = (item *) malloc (nbItem * sizeof(item));
    for(i = 0;i < nbItem;i++)
    {
        prob[i].p = probp[i];
        prob[i].w = probw[i];
        prob[i].i = probi[i];
        prob[i].x = probx[i];
    }
    return combo(prob,&(prob[nbItem-1]),capa,0,maxZ,1,1);
}
