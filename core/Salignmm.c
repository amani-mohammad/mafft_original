//S alignment methods
#include "mltaln.h"
#include "dp.h"

#define MACHIGAI 0
#define OUTGAP0TRY 0
#define DEBUG 0
#define XXXXXXX    0
#define USE_PENALTY_EX  1
#define FASTMATCHCALC 1
#define SLOW 0


static TLS double **impmtx = NULL;
static TLS int impalloclen = 0;
double imp_match_out_sc( int i1, int j1 )
{
//	fprintf( stderr, "imp+match = %f\n", impmtx[i1][j1] * fastathreshold );
//	fprintf( stderr, "val = %f\n", impmtx[i1][j1] );
	return( impmtx[i1][j1] );
}

#if 0
static void imp_match_out_vead_gapmap( double *imp, int i1, int lgth2, int *gapmap2 )
{
#if FASTMATCHCALC
	double *pt = impmtx[i1];
	int *gapmappt = gapmap2;
	while( lgth2-- )
		*imp++ += pt[*gapmappt++];
#else
	int j;
	double *pt = impmtx[i1];
	for( j=0; j<lgth2; j++ )
		*imp++ += pt[gapmap2[j]];
#endif
}


static void imp_match_out_vead_tate_gapmap( double *imp, int j1, int lgth1, int *gapmap1 )
{
#if FASTMATCHCALC
	int *gapmappt = gapmap1;
	while( lgth1-- )
		*imp++ += impmtx[*gapmappt++][j1];
#else
	int i;
	for( i=0; i<lgth1; i++ )
		*imp++ += impmtx[gapmap1[i]][j1];
#endif
}
#endif

static void imp_match_out_vead( double *imp, int i1, int lgth2 )
{
#if FASTMATCHCALC 
	double *pt = impmtx[i1];
	while( lgth2-- )
		*imp++ += *pt++;
#else
	int j;
	double *pt = impmtx[i1];
	for( j=0; j<lgth2; j++ )
		*imp++ += pt[j];
#endif
}
static void imp_match_out_vead_tate( double *imp, int j1, int lgth1 )
{
	int i;
	for( i=0; i<lgth1; i++ )
		*imp++ += impmtx[i][j1];
}

void imp_rna( int nseq1, int nseq2, char **seq1, char **seq2, double *eff1, double *eff2, RNApair ***grouprna1, RNApair ***grouprna2, int *gapmap1, int *gapmap2, RNApair *pair )
{
	foldrna( nseq1, nseq2, seq1, seq2, eff1, eff2, grouprna1, grouprna2, impmtx, gapmap1, gapmap2, pair );
}

void imp_match_init_strict( double *imp, int clus1, int clus2, int lgth1, int lgth2, char **seq1, char **seq2, double *eff1, double *eff2, double *eff1_kozo, double *eff2_kozo, LocalHom ***localhom, char *swaplist, int forscore, int *orinum1, int *orinum2 )
{
//	int i, j, k1, k2, tmpint, start1, start2, end1, end2;
//	double effij;
//	double effij_kozo;
//	double effijx;
//	char *pt, *pt1, *pt2;
//	static TLS char *nocount1 = NULL;
//	static TLS char *nocount2 = NULL;
//	LocalHom *tmpptr;

	if( seq1 == NULL )
	{
		if( impmtx ) FreeFloatMtx( impmtx );
		impmtx = NULL;
//		if( nocount1 ) free( nocount1 );
//		nocount1 = NULL;
//		if( nocount2 ) free( nocount2 );
//		nocount2 = NULL;
		
		return;
	}

	if( impalloclen < lgth1 + 2 || impalloclen < lgth2 + 2 )
	{
		if( impmtx ) FreeFloatMtx( impmtx );
//		if( nocount1 ) free( nocount1 );
//		if( nocount2 ) free( nocount2 );
		impalloclen = MAX( lgth1, lgth2 ) + 2;
		impmtx = AllocateFloatMtx( impalloclen, impalloclen );
//		nocount1 = AllocateCharVec( impalloclen );
//		nocount2 = AllocateCharVec( impalloclen );
	}

	fillimp( impmtx, imp, clus1, clus2, lgth1, lgth2, seq1, seq2, eff1, eff2, eff1_kozo, eff2_kozo, localhom, swaplist, forscore, orinum1, orinum2 );
}


static void match_calc_del( int **which, double ***matrices, double *match, int n1, char **seq1, double *eff1, int n2, char **seq2, double *eff2, int i1, int lgth2, int mid, int nmask, int *mask1, int *mask2 ) 
{
// osoi!
	int i, j, k, m;
	int c1, c2;
//	fprintf( stderr, "\nmatch_calc_dynamicmtx... %d", i1 );
//	fprintf( stderr, "\nseq1[0]=%s\n", seq1[0] );
//	fprintf( stderr, "\nseq2[0]=%s\n", seq2[0] );
//	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
//	{
//		if( flip ) reporterr( "in match_calc_slow, which[%d][%d] = %d\n", j, i, which[j][i] );
//		else       reporterr( "in match_calc_slow, which[%d][%d] = %d\n", i, j, which[i][j] );
//	}
	for( k=0; k<lgth2; k++ )
	{
		for( m=0; m<nmask; m++ )
		{
			i = mask1[m];
			j = mask2[m];
//			reporterr( "Deleting %d-%d (c=%d)\n", i, j, mid );
//			if( k==0 ) fprintf( stderr, "pairoffset[%d][%d] = %f\n", i, j, po );
			c1 = amino_n[(unsigned int)seq1[i][i1]];
			c2 = amino_n[(unsigned int)seq2[j][k]];
//			reporterr( "k=%d, c1=%d, c2=%d, seq1[i][i1]=%c, seq2[%d][%d]=%c\n", k, c1, c2, seq1[i][i1], j, k, seq2[j][k] );
			if( seq1[i][i1] == '-' || seq2[j][k] == '-' ) continue;
			if( c1 < 0 || c2 < 0 ) continue;
//			fprintf( stderr, "c1=%d, c2=%d\n", c1, c2 );
//			fprintf( stderr, "match[k] = %f -> ", match[k], mid );
			match[k] -= matrices[mid][c1][c2] * eff1[i] * eff2[j];
//			fprintf( stderr, "match[k] = %f (mid=%d)\n", match[k], mid );
		}
	}
//	fprintf( stderr, "done\n" );
	return;
}


#if SLOW
static void match_calc_slow( int **which, double ***matrices, double *match, int n1, char **seq1, double *eff1, int n2, char **seq2, double *eff2, int i1, int lgth2, double **doublework, int **intwork, int initialize, int flip ) 
{
// osoi!
	int i, j, k;
	int c1, c2;
	int mid;
//	fprintf( stderr, "\nmatch_calc_dynamicmtx... %d", i1 );
//	fprintf( stderr, "\nseq1[0]=%s\n", seq1[0] );
//	fprintf( stderr, "\nseq2[0]=%s\n", seq2[0] );
//	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
//	{
//		if( flip ) reporterr( "in match_calc_slow, which[%d][%d] = %d\n", j, i, which[j][i] );
//		else       reporterr( "in match_calc_slow, which[%d][%d] = %d\n", i, j, which[i][j] );
//	}
	for( k=0; k<lgth2; k++ )
	{
		match[k] = 0.0;
		for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
		{
			if( flip ) mid = which[j][i];
			else       mid = which[i][j];
//			if( k==0 ) fprintf( stderr, "pairoffset[%d][%d] = %f\n", i, j, po );
			c1 = amino_n[(unsigned int)seq1[i][i1]];
			c2 = amino_n[(unsigned int)seq2[j][k]];
			if( seq1[i][i1] == '-' || seq2[j][k] == '-' ) continue;
			if( c1 < 0 || c2 < 0 ) continue;
//			fprintf( stderr, "c1=%d, c2=%d\n", c1, c2 );
			if( flip ) 
				match[k] += matrices[mid][c1][c2] * eff1[i] * eff2[j];
			else
				match[k] += matrices[mid][c1][c2] * eff1[i] * eff2[j];
//			fprintf( stderr, "match[k] = %f (which=%d)\n", match[k], mid );
		}
	}
//	fprintf( stderr, "done\n" );
	return;
}
#endif

static void fillzero( double *s, int l )
{
	while( l-- ) *s++ = 0.0;
}

static void match_calc_add( double **scoreingmtx, double *match, double **cpmx1, double **cpmx2, int i1, int lgth2, double **doublework, int **intwork, int initialize )
{
#if FASTMATCHCALC
//	fprintf( stderr, "\nmatch_calc... %d", i1 );
	int j, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *matchpt, *cpmxpdpt, **cpmxpdptpt;
	int *cpmxpdnpt, **cpmxpdnptpt;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
	if( initialize )
	{
		int count = 0;
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[j][count] = cpmx2[l][j];
					cpmxpdn[j][count] = l;
					count++;
				}
			}
			cpmxpdn[j][count] = -1;
		}
	}

	{
		for( l=0; l<nalphabets; l++ )
		{
			scarr[l] = 0.0;
			for( j=0; j<nalphabets; j++ )
//				scarr[l] += n_dis[j][l] * cpmx1[j][i1];
//				scarr[l] += n_dis_consweight_multi[j][l] * cpmx1[j][i1];
				scarr[l] += scoreingmtx[j][l] * cpmx1[j][i1];
		}
		matchpt = match;
		cpmxpdnptpt = cpmxpdn;
		cpmxpdptpt = cpmxpd;
		while( lgth2-- )
		{
//			*matchpt = 0.0;
			cpmxpdnpt = *cpmxpdnptpt++;
			cpmxpdpt = *cpmxpdptpt++;
			while( *cpmxpdnpt>-1 )
				*matchpt += scarr[*cpmxpdnpt++] * *cpmxpdpt++;
			matchpt++;
		} 
	}
	free( scarr );
//	fprintf( stderr, "done\n" );
#else
	int j, k, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
// simple
	if( initialize )
	{
		int count = 0;
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[count][j] = cpmx2[l][j];
					cpmxpdn[count][j] = l;
					count++;
				}
			}
			cpmxpdn[count][j] = -1;
		}
	}
	for( l=0; l<nalphabets; l++ )
	{
		scarr[l] = 0.0;
		for( k=0; k<nalphabets; k++ )
//			scarr[l] += n_dis[k][l] * cpmx1[k][i1];
//			scarr[l] += n_dis_consweight_multi[k][l] * cpmx1[k][i1];
			scarr[l] += scoreingmtx[k][l] * cpmx1[k][i1];
	}
	for( j=0; j<lgth2; j++ )
	{
		match[j] = 0.0;
		for( k=0; cpmxpdn[k][j]>-1; k++ )
			match[j] += scarr[cpmxpdn[k][j]] * cpmxpd[k][j];
	} 
	free( scarr );
#endif
}
static void match_calc( double **n_dynamicmtx, double *match, double **cpmx1, double **cpmx2, int i1, int lgth2, double **doublework, int **intwork, int initialize )
{
#if FASTMATCHCALC
//	fprintf( stderr, "\nmatch_calc... %d", i1 );
	int j, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *matchpt, *cpmxpdpt, **cpmxpdptpt;
	int *cpmxpdnpt, **cpmxpdnptpt;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
	if( initialize )
	{
		int count = 0;
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[j][count] = cpmx2[l][j];
					cpmxpdn[j][count] = l;
					count++;
				}
			}
			cpmxpdn[j][count] = -1;
		}
	}

	{
		for( l=0; l<nalphabets; l++ )
		{
			scarr[l] = 0.0;
			for( j=0; j<nalphabets; j++ )
//				scarr[l] += n_dis[j][l] * cpmx1[j][i1];
//				scarr[l] += n_dis_consweight_multi[j][l] * cpmx1[j][i1];
				scarr[l] += n_dynamicmtx[j][l] * cpmx1[j][i1];
		}
		matchpt = match;
		cpmxpdnptpt = cpmxpdn;
		cpmxpdptpt = cpmxpd;
		while( lgth2-- )
		{
			*matchpt = 0.0;
			cpmxpdnpt = *cpmxpdnptpt++;
			cpmxpdpt = *cpmxpdptpt++;
			while( *cpmxpdnpt>-1 )
				*matchpt += scarr[*cpmxpdnpt++] * *cpmxpdpt++;
			matchpt++;
		} 
	}
	free( scarr );
//	fprintf( stderr, "done\n" );
#else
	int j, k, l;
//	double scarr[26];
	double **cpmxpd = doublework;
	int **cpmxpdn = intwork;
	double *scarr;
	scarr = calloc( nalphabets, sizeof( double ) );
// simple
	if( initialize )
	{
		int count = 0;
		for( j=0; j<lgth2; j++ )
		{
			count = 0;
			for( l=0; l<nalphabets; l++ )
			{
				if( cpmx2[l][j] )
				{
					cpmxpd[count][j] = cpmx2[l][j];
					cpmxpdn[count][j] = l;
					count++;
				}
			}
			cpmxpdn[count][j] = -1;
		}
	}
	for( l=0; l<nalphabets; l++ )
	{
		scarr[l] = 0.0;
		for( k=0; k<nalphabets; k++ )
//			scarr[l] += n_dis[k][l] * cpmx1[k][i1];
//			scarr[l] += n_dis_consweight_multi[k][l] * cpmx1[k][i1];
			scarr[l] += n_dynamicmtx[k][l] * cpmx1[k][i1];
	}
	for( j=0; j<lgth2; j++ )
	{
		match[j] = 0.0;
		for( k=0; cpmxpdn[k][j]>-1; k++ )
			match[j] += scarr[cpmxpdn[k][j]] * cpmxpd[k][j];
	} 
	free( scarr );
#endif
}

static void Atracking_localhom( double *impwmpt, double *lasthorizontalw, double *lastverticalw, 
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        int **ijp, int icyc, int jcyc,
						int *warpis, int *warpjs, int warpbase,
						int *ngap1, int *ngap2, int reuseprofiles )
{
	int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, limk;
	double wm;
	char *gaptable1, *gt1bk;
	char *gaptable2, *gt2bk;
	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );
	gt1bk = AllocateCharVec( lgth1+lgth2+1 );
	gt2bk = AllocateCharVec( lgth1+lgth2+1 );

#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif
 
	if( outgap == 1 )
		;
	else
	{
		wm = lastverticalw[0];
		for( i=0; i<lgth1; i++ )
		{
			if( lastverticalw[i] >= wm )
			{
				wm = lastverticalw[i];
				iin = i; jin = lgth2-1;
				ijp[lgth1][lgth2] = +( lgth1 - i );
			}
		}
		for( j=0; j<lgth2; j++ )
		{
			if( lasthorizontalw[j] >= wm )
			{
				wm = lasthorizontalw[j];
				iin = lgth1-1; jin = j;
				ijp[lgth1][lgth2] = -( lgth2 - j );
			}
		}
	}

    for( i=0; i<lgth1+1; i++ ) 
    {
        ijp[i][0] = i + 1;
    }
    for( j=0; j<lgth2+1; j++ ) 
    {
        ijp[0][j] = -( j + 1 );
    }

	gaptable1 = gt1bk + lgth1+lgth2;
	*gaptable1 = 0;
	gaptable2 = gt2bk + lgth1+lgth2;
	*gaptable2 = 0;
	*ngap1 = *ngap2 = 0;

	iin = lgth1; jin = lgth2;
	limk = lgth1+lgth2 + 1;
	*impwmpt = 0.0;
	for( k=0; k<limk; k++ ) 
	{
		if( ijp[iin][jin] >= warpbase )
		{
			ifi = warpis[ijp[iin][jin]-warpbase];
			jfi = warpjs[ijp[iin][jin]-warpbase];
		}
		else if( ijp[iin][jin] < 0 ) 
		{
			ifi = iin-1; jfi = jin+ijp[iin][jin];
		}
		else if( ijp[iin][jin] > 0 )
		{
			ifi = iin-ijp[iin][jin]; jfi = jin-1;
		}
		else
		{
			ifi = iin-1; jfi = jin-1;
		}
		if( ifi == -warpbase && jfi == -warpbase )
		{
			l = iin;
			while( --l >= 0 )
			{
				*--gaptable1 = 'o';
				*--gaptable2 = '-';
				k++;
				*ngap2 = 1;
			}
			l= jin;
			while( --l >= 0 )
			{
				*--gaptable1 = '-';
				*--gaptable2 = 'o';
				k++;
				*ngap1 = 1;
			}
			break;
		}
		else
		{
			l = iin - ifi;
			while( --l ) 
			{
				*--gaptable1 = 'o';
				*--gaptable2 = '-';
				k++;
				*ngap2 = 1;
			}
			l= jin - jfi;
			while( --l )
			{
				*--gaptable1 = '-';
				*--gaptable2 = 'o';
				k++;
				*ngap1 = 1;
			}
		}
		if( iin == lgth1 || jin == lgth2 )
			;
		else
		{
			*impwmpt += (double)imp_match_out_sc( iin, jin );

//		fprintf( stderr, "impwm = %f (iin=%d, jin=%d) seq1=%c, seq2=%c\n", *impwmpt, iin, jin, seq1[0][iin], seq2[0][jin] );
		}
		if( iin <= 0 || jin <= 0 ) break;
		*--gaptable1 = 'o';
		*--gaptable2 = 'o';
		k++;
		iin = ifi; jin = jfi;
	}

	if( *ngap1 || !reuseprofiles )
		for( i=0; i<icyc; i++ ) gapireru( mseq1[i], seq1[i], gaptable1 );

	if( *ngap2 || !reuseprofiles )
		for( j=0; j<jcyc; j++ ) gapireru( mseq2[j], seq2[j], gaptable2 );

	free( gt1bk );
	free( gt2bk );
}

static double Atracking( double *lasthorizontalw, double *lastverticalw, 
						char **seq1, char **seq2, 
                        char **mseq1, char **mseq2, 
                        int **ijp, int icyc, int jcyc,
						int tailgp,
						int *warpis, int *warpjs, int warpbase,
						int *ngap1, int *ngap2, int reuseprofiles )
{
	int i, j, l, iin, jin, ifi, jfi, lgth1, lgth2, k, limk;
	double wm;
	char *gaptable1, *gt1bk;
	char *gaptable2, *gt2bk;
	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );

	gt1bk = AllocateCharVec( lgth1+lgth2+1 );
	gt2bk = AllocateCharVec( lgth1+lgth2+1 );

#if 0
	for( i=0; i<lgth1; i++ ) 
	{
		fprintf( stderr, "lastverticalw[%d] = %f\n", i, lastverticalw[i] );
	}
#endif
 
	if( tailgp == 1 )
		;
	else
	{
		wm = lastverticalw[0];
		for( i=0; i<lgth1; i++ )
		{
			if( lastverticalw[i] >= wm )
			{
				wm = lastverticalw[i];
				iin = i; jin = lgth2-1;
				ijp[lgth1][lgth2] = +( lgth1 - i );
			}
		}
		for( j=0; j<lgth2; j++ )
		{
			if( lasthorizontalw[j] >= wm )
			{
				wm = lasthorizontalw[j];
				iin = lgth1-1; jin = j;
				ijp[lgth1][lgth2] = -( lgth2 - j );
			}
		}
	}

    for( i=0; i<lgth1+1; i++ ) 
    {
        ijp[i][0] = i + 1;
    }
    for( j=0; j<lgth2+1; j++ ) 
    {
        ijp[0][j] = -( j + 1 );
    }

	gaptable1 = gt1bk + lgth1+lgth2;
	*gaptable1 = 0;
	gaptable2 = gt2bk + lgth1+lgth2;
	*gaptable2 = 0;

	*ngap1 = *ngap2 = 0;
	iin = lgth1; jin = lgth2;
	limk = lgth1+lgth2 + 1;
	for( k=0; k<limk; k++ ) 
	{
		if( ijp[iin][jin] >= warpbase )
		{
			ifi = warpis[ijp[iin][jin]-warpbase];
			jfi = warpjs[ijp[iin][jin]-warpbase];
		}
		else if( ijp[iin][jin] < 0 ) 
		{
			ifi = iin-1; jfi = jin+ijp[iin][jin];
		}
		else if( ijp[iin][jin] > 0 )
		{
			ifi = iin-ijp[iin][jin]; jfi = jin-1;
		}
		else
		{
			ifi = iin-1; jfi = jin-1;
		}

		if( ifi == -warpbase && jfi == -warpbase )
		{
			l = iin;
			while( --l >= 0 )
			{
				*--gaptable1 = 'o';
				*--gaptable2 = '-';
				k++;
				*ngap2=1;
			}
			l= jin;
			while( --l >= 0 )
			{
				*--gaptable1 = '-';
				*--gaptable2 = 'o';
				k++;
				*ngap1=1;
			}
			break;
		}
		else
		{
			l = iin - ifi;
			while( --l ) 
			{
				*--gaptable1 = 'o';
				*--gaptable2 = '-';
				k++;
				*ngap2=1;
			}
			l= jin - jfi;
			while( --l )
			{
				*--gaptable1 = '-';
				*--gaptable2 = 'o';
				k++;
				*ngap1=1;
			}
		}
		if( iin <= 0 || jin <= 0 ) break;
		*--gaptable1 = 'o';
		*--gaptable2 = 'o';
		k++;
		iin = ifi; jin = jfi;
	}

//	reporterr( "gaptable1=%s\n", gaptable1 );
//	reporterr( "gaptable2=%s\n", gaptable2 );


	if( *ngap1 || !reuseprofiles )
		for( i=0; i<icyc; i++ ) gapireru( mseq1[i], seq1[i], gaptable1 );

	if( *ngap2 || !reuseprofiles )
		for( j=0; j<jcyc; j++ ) gapireru( mseq2[j], seq2[j], gaptable2 );

	free( gt1bk );
	free( gt2bk );

	return( 0.0 );
}

double A__align( double **n_dynamicmtx, char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen, LocalHom ***localhom, double *impmatch, char *sgap1, char *sgap2, char *egap1, char *egap2, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp, int firstmem, int calledbyfulltreebase )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
	int reuseprofiles;
	static TLS int previousfirstlen; // 2016/Feb/1 // MEMBER NO CHECK GA HITSUYOU!!!!
	static TLS int previousicyc; // 2016/Feb/1 // MEMBER NO CHECK GA HITSUYOU!!!!
	static TLS int previousfirstmem;
	static TLS int previouscall;
	int ngap1, ngap2;
//	int k;
	register int i, j;
	int lasti, lastj;      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lgth1, lgth2;
	int resultlen;
	double wm = 0.0;   /* int ?????? */
	double g;
	double *currentw, *previousw;
//	double fpenalty = (double)penalty;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex;
#endif
#if 1
	double *wtmp;
	int *ijppt;
	double *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	static TLS double mi, *m;
	static TLS int **ijp;
	static TLS int mpi, *mp;
	static TLS double *w1, *w2;
	static TLS double *match;
	static TLS double *initverticalw;    /* kufuu sureba iranai */
	static TLS double *lastverticalw;    /* kufuu sureba iranai */
	static TLS char **mseq1;
	static TLS char **mseq2;
	static TLS char **mseq;
	static TLS double *ogcp1, *ogcp1o;
	static TLS double *ogcp2;
	static TLS double *fgcp1, *fgcp1o;
	static TLS double *fgcp2;
	static TLS double **cpmx1;
	static TLS double **cpmx2;
	static TLS int **intwork;
	static TLS double **doublework;
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double *gapfreq1;
	static TLS double *gapfreq2;
	double fpenalty = (double)penalty;
	double fpenalty_shift = (double)penalty_shift;
	double *fgcp2pt;
	double *ogcp2pt;
	double fgcp1va;
	double ogcp1va;
	double *gf2pt;
	double *gf2ptpre;
	double gf1va;
	double gf1vapre;
	double headgapfreq1;
	double headgapfreq2;

	int *warpis = NULL;
	int *warpjs = NULL;
	int *warpi = NULL;
	int *warpj = NULL;
	int *prevwarpi = NULL;
	int *prevwarpj = NULL;
	double *wmrecords = NULL;
	double *prevwmrecords = NULL;
	int warpn = 0;
	int warpbase;
	double curm = 0.0;
	double *wmrecordspt, *wmrecords1pt, *prevwmrecordspt;
	int *warpipt, *warpjpt;

//	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s, %f\n", seq1[i], eff1[i] );
//	for( i=0; i<jcyc; i++ ) fprintf( stderr, "%s, %f\n", seq2[i], eff2[i] );


	if( seq1 == NULL )
	{
		if( orlgth1 )
		{
//			fprintf( stderr, "## Freeing local arrays in A__align\n" );
			orlgth1 = 0;
			orlgth2 = 0;

			imp_match_init_strict( NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, NULL );

			free( mseq1 );
			free( mseq2 );
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );

			FreeFloatVec( ogcp1 );
			FreeFloatVec( ogcp1o );
			FreeFloatVec( ogcp2 );
			FreeFloatVec( fgcp1 );
			FreeFloatVec( fgcp1o );
			FreeFloatVec( fgcp2 );


			FreeFloatMtx( cpmx1 );
			FreeFloatMtx( cpmx2 );

			FreeFloatVec( gapfreq1 );
			FreeFloatVec( gapfreq2 );

			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );

		}
		else
		{
//			fprintf( stderr, "## Not allocated\n" );
		}
		return( 0.0 );
	}


	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );
#if 0
	if( lgth1 == 0 || lgth2 == 0 )
	{
		fprintf( stderr, "WARNING (Aalignmm): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif
	if( lgth1 == 0 && lgth2 == 0 )
		return( 0.0 );

	if( lgth1 == 0 )
	{
		for( i=0; i<icyc; i++ )
		{
			j = lgth2;
			seq1[i][j] = 0;
			while( j ) seq1[i][--j] = *newgapstr;
//			fprintf( stderr, "seq1[i] = %s\n", seq1[i] );
		}
		return( 0.0 );
	}

	if( lgth2 == 0 )
	{
		for( i=0; i<jcyc; i++ )
		{
			j = lgth1;
			seq2[i][j] = 0;
			while( j ) seq2[i][--j] = *newgapstr;
//			fprintf( stderr, "seq2[i] = %s\n", seq2[i] );
		}
		return( 0.0 );
	}

	warpbase = lgth1 + lgth2;
	warpis = NULL;
	warpjs = NULL;
	warpn = 0;



	if( trywarp )
	{
//		fprintf( stderr, "IN A__align, penalty_shift = %d\n", penalty_shift );
		if( headgp == 0 || tailgp == 0 )
		{
			fprintf( stderr, "At present, headgp and tailgp must be 1 to allow shift.\n" );
			exit( 1 );
		}
		wmrecords = AllocateFloatVec( lgth2+1 );
		warpi = AllocateIntVec( lgth2+1 );
		warpj = AllocateIntVec( lgth2+1 );
		prevwmrecords = AllocateFloatVec( lgth2+1 );
		prevwarpi = AllocateIntVec( lgth2+1 );
		prevwarpj = AllocateIntVec( lgth2+1 );
		for( i=0; i<lgth2+1; i++ ) wmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwarpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) prevwarpj[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpj[i] = -warpbase;
	}


#if 0
	fprintf( stderr, "####  eff in SA+++align\n" );
	fprintf( stderr, "####  seq1[0] = %s\n", seq1[0] );
	fprintf( stderr, "####  strlen( seq1[0] ) = %d\n", strlen( seq1[0] ) );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "eff1[%d] = %f\n", i, eff1[i] );
	fprintf( stderr, "####  seq2[0] = %s\n", seq2[0] );
	fprintf( stderr, "####  strlen( seq2[0] ) = %d\n", strlen( seq2[0] ) );
	for( i=0; i<jcyc; i++ ) fprintf( stderr, "eff2[%d] = %f\n", i, eff2[i] );
#endif
	if( orlgth1 == 0 )
	{
		mseq1 = AllocateCharMtx( njob, 0 );
		mseq2 = AllocateCharMtx( njob, 0 );
	}



	if( lgth1 > orlgth1 || lgth2 > orlgth2 )
	{
		int ll1, ll2;


		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );

			FreeFloatVec( ogcp1 );
			FreeFloatVec( ogcp1o );
			FreeFloatVec( ogcp2 );
			FreeFloatVec( fgcp1 );
			FreeFloatVec( fgcp1o );
			FreeFloatVec( fgcp2 );


			FreeFloatMtx( cpmx1 );
			FreeFloatMtx( cpmx2 );

			FreeFloatVec( gapfreq1 );
			FreeFloatVec( gapfreq2 );

			FreeFloatMtx( doublework );
			FreeIntMtx( intwork );
		}

		ll1 = MAX( (int)(1.3*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.3*lgth2), orlgth2 ) + 100;

#if DEBUG
		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );
#endif

		w1 = AllocateFloatVec( ll2+2 );
		w2 = AllocateFloatVec( ll2+2 );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = AllocateFloatVec( ll2+2 );
		mp = AllocateIntVec( ll2+2 );

		mseq = AllocateCharMtx( njob, ll1+ll2 );

		ogcp1 = AllocateFloatVec( ll1+2 );
		ogcp1o = AllocateFloatVec( ll1+2 );
		ogcp2 = AllocateFloatVec( ll2+2 );
		fgcp1 = AllocateFloatVec( ll1+2 );
		fgcp1o = AllocateFloatVec( ll1+2 );
		fgcp2 = AllocateFloatVec( ll2+2 );

		cpmx1 = AllocateFloatMtx( nalphabets, ll1+2 );
		cpmx2 = AllocateFloatMtx( nalphabets, ll2+2 );
		previousfirstlen = -1;
		previousicyc = -1;

		gapfreq1 = AllocateFloatVec( ll1+2 );
		gapfreq2 = AllocateFloatVec( ll2+2 );

#if FASTMATCHCALC
		doublework = AllocateFloatMtx( MAX( ll1, ll2 )+2, nalphabets ); 
		intwork = AllocateIntMtx( MAX( ll1, ll2 )+2, nalphabets+1 ); 
#else
		doublework = AllocateFloatMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
		intwork = AllocateIntMtx( nalphabets, MAX( ll1, ll2 )+2 ); 
#endif

#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}


	for( i=0; i<icyc; i++ )
	{
		mseq1[i] = mseq[i];
		seq1[i][lgth1] = 0;
	}
	for( j=0; j<jcyc; j++ )
	{
		mseq2[j] = mseq[icyc+j];
		seq2[j][lgth2] = 0;
	}


	if( orlgth1 > commonAlloc1 || orlgth2 > commonAlloc2 )
	{
		int ll1, ll2;

		if( commonAlloc1 && commonAlloc2 )
		{
			FreeIntMtx( commonIP );
		}

		ll1 = MAX( orlgth1, commonAlloc1 );
		ll2 = MAX( orlgth2, commonAlloc2 );

#if DEBUG
		fprintf( stderr, "\n\ntrying to allocate %dx%d matrices ... ", ll1+1, ll2+1 );
#endif

		commonIP = AllocateIntMtx( ll1+10, ll2+10 );

#if DEBUG
		fprintf( stderr, "succeeded\n\n" );
#endif

		commonAlloc1 = ll1;
		commonAlloc2 = ll2;
	}
	ijp = commonIP;

#if 0
	{
		double t = 0.0;
		for( i=0; i<icyc; i++ )
			t += eff1[i];
	fprintf( stderr, "## totaleff = %f\n", t );
	}
#endif

	if( calledbyfulltreebase && previouscall && firstmem >= 0 && firstmem == previousfirstmem && lgth1 == previousfirstlen && nthread == 0 && icyc == previousicyc+1 ) // mouhitotsu jouken ga iru. firstmem >= 0 && firstmem == previousfirstmem && calledbytreebase && prevcalledbytreebase
		reuseprofiles = 1;
	else
		reuseprofiles = 0;


	if( reuseprofiles )
	{
//		reporterr( "r" );
		cpmx_calc_add( seq1, cpmx1, eff1, lgth1, icyc );
	}
	else
	{
//		reporterr( "n" );
		cpmx_calc_new( seq1, cpmx1, eff1, lgth1, icyc );
	}
	cpmx_calc_new( seq2, cpmx2, eff2, lgth2, jcyc );

	if( sgap1 )
	{
		new_OpeningGapCount( ogcp1o, icyc, seq1, eff1, lgth1, sgap1 );
		new_FinalGapCount( fgcp1o, icyc, seq1, eff1, lgth1, egap1 );

		new_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2, sgap2 );
		new_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2, egap2 );

		outgapcount( &headgapfreq1, icyc, sgap1, eff1 );
		outgapcount( &headgapfreq2, jcyc, sgap2, eff2 );
		outgapcount( gapfreq1+lgth1, icyc, egap1, eff1 );
		outgapcount( gapfreq2+lgth2, jcyc, egap2, eff2 );
	}
	else
	{
		if( reuseprofiles )
		{
			st_OpeningGapAdd( ogcp1o, icyc, seq1, eff1, lgth1 );
			st_FinalGapAdd( fgcp1o, icyc, seq1, eff1, lgth1 );
		}
		else
		{
			st_OpeningGapCount( ogcp1o, icyc, seq1, eff1, lgth1 );
			st_FinalGapCount( fgcp1o, icyc, seq1, eff1, lgth1 );
		}

		st_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2 );
		st_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2 );

		headgapfreq1 = 0.0;
		headgapfreq2 = 0.0;
		gapfreq1[lgth1] = 0.0;
		gapfreq2[lgth2] = 0.0;
	}

	if( legacygapcost == 0 )
	{
		if( reuseprofiles )
			gapcountadd( gapfreq1, seq1, icyc, eff1, lgth1 );
		else
			gapcountf( gapfreq1, seq1, icyc, eff1, lgth1 );

		gapcountf( gapfreq2, seq2, jcyc, eff2, lgth2 );

		for( i=0; i<lgth1+1; i++ ) gapfreq1[i] = 1.0 - gapfreq1[i];
		for( i=0; i<lgth2+1; i++ ) gapfreq2[i] = 1.0 - gapfreq2[i];
		headgapfreq1 = 1.0 - headgapfreq1;
		headgapfreq2 = 1.0 - headgapfreq2;
	}
	else
	{
		for( i=0; i<lgth1+1; i++ ) gapfreq1[i] = 1.0;
		for( i=0; i<lgth2+1; i++ ) gapfreq2[i] = 1.0;
		headgapfreq1 = 1.0;
		headgapfreq2 = 1.0;
	}

#if 0
	fprintf( stderr, "\ngapfreq1[] =" );
	for( i=0; i<lgth1; i++ ) fprintf( stderr, "%5.2f ", gapfreq1[i] );
	fprintf( stderr, "\n" );

	fprintf( stderr, "\ngapfreq2[] =" );
	for( i=0; i<lgth2; i++ ) fprintf( stderr, "%5.2f ", gapfreq2[i] );
	fprintf( stderr, "\n" );
#endif
	

	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 0.5 * ( 1.0 - ogcp1o[i] ) * fpenalty * ( gapfreq1[i] );
		fgcp1[i] = 0.5 * ( 1.0 - fgcp1o[i] ) * fpenalty * ( gapfreq1[i] );
	}

	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 0.5 * ( 1.0 - ogcp2[i] ) * fpenalty * ( gapfreq2[i] );
		fgcp2[i] = 0.5 * ( 1.0 - fgcp2[i] ) * fpenalty * ( gapfreq2[i] );
	}
#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;

	match_calc( n_dynamicmtx, initverticalw, cpmx2, cpmx1, 0, lgth1, doublework, intwork, 1 );
	if( localhom )
		imp_match_out_vead_tate( initverticalw, 0, lgth1 ); // 060306

	match_calc( n_dynamicmtx, currentw, cpmx1, cpmx2, 0, lgth2, doublework, intwork, 1 );
	if( localhom )
		imp_match_out_vead( currentw, 0, lgth2 ); // 060306
#if 0 // -> tbfast.c
	if( localhom )
		imp_match_calc( n_dynamicmtx, currentw, icyc, jcyc, lgth1, lgth2, seq1, seq2, eff1, eff2, localhom, 1, 0 );

#endif

	if( headgp == 1 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
//			initverticalw[i] += ( ogcp1[0] + fgcp1[i-1] ) ;
			initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2[0] ) ;
		}
		for( j=1; j<lgth2+1; j++ )
		{
//			currentw[j] += ( ogcp2[0] + fgcp2[j-1] ) ;
			currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1[0] ) ;
		}
	}
#if OUTGAP0TRY
	else
	{
		fprintf( stderr, "offset = %d\n", offset );
		for( j=1; j<lgth2+1; j++ )
			currentw[j] -= offset * j / 2.0;
		for( i=1; i<lgth1+1; i++ )
			initverticalw[i] -= offset * i / 2.0;
	}
#endif
#if 0
	fprintf( stderr, "\n   " );
	for( j=0; j<lgth2+1; j++ ) fprintf( stderr, "    %c ", seq2[0][j] );
	fprintf( stderr, "\n%c ", seq1[0][0]  );
	for( j=0; j<lgth2+1; j++ )
	{
		fprintf( stderr, "%5.0f ", currentw[j] );
	}
	fprintf( stderr, "\n"  );
#endif


	for( j=1; j<lgth2+1; ++j ) 
	{
//		m[j] = currentw[j-1] + ogcp1[1]; mp[j] = 0;
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2[j-1]; mp[j] = 0;;
	}
	if( lgth2 == 0 )
		lastverticalw[0] = 0.0; // Falign kara yobaretatoki kounarukanousei ari
	else
		lastverticalw[0] = currentw[lgth2-1];

	if( tailgp ) lasti = lgth1+1; else lasti = lgth1;
	lastj = lgth2+1;

#if XXXXXXX
fprintf( stderr, "currentw = \n" );
for( i=0; i<lgth1+1; i++ )
{
	fprintf( stderr, "%5.2f ", currentw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "initverticalw = \n" );
for( i=0; i<lgth2+1; i++ )
{
	fprintf( stderr, "%5.2f ", initverticalw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "fcgp\n" );
for( i=0; i<lgth1; i++ ) 
	fprintf( stderr, "fgcp1[%d]=%f\n", i, ogcp1[i] );
for( i=0; i<lgth2; i++ ) 
	fprintf( stderr, "fgcp2[%d]=%f\n", i, ogcp2[i] );
#endif

	for( i=1; i<lasti; i++ )
	{

#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
			*chudanres = 1;
			return( -1.0 );
		}
#endif
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

		match_calc( n_dynamicmtx, currentw, cpmx1, cpmx2, i, lgth2, doublework, intwork, 0 );
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
		if( localhom )
		{
//			fprintf( stderr, "Calling imp_match_calc (o) lgth = %d, i = %d\n", lgth1, i );
#if  0
			imp_match_out_vead( currentw, i, lgth2 );
#else
			imp_match_out_vead( currentw, i, lgth2 );
#endif
		}
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
		currentw[0] = initverticalw[i];

#if 0
		fprintf( stderr, "%c ", seq1[0][i] );
		for( j=0; j<lgth2+1; j++ )
		{
			fprintf( stderr, "%5.0f ", currentw[j] );
		}
		fprintf( stderr, "\n"  );
#endif
	
//		mi = previousw[0] + ogcp2[1]; mpi = 0;
		mi = previousw[0] + ogcp2[1] * gapfreq1[i-1]; mpi=0;
		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
		fgcp2pt = fgcp2;
		ogcp2pt = ogcp2 + 1;
		fgcp1va = fgcp1[i-1];
		ogcp1va = ogcp1[i];
		gf1va = gapfreq1[i];
		gf1vapre = gapfreq1[i-1];
		gf2pt = gapfreq2+1;
		gf2ptpre = gapfreq2;

		if( trywarp )
		{
			prevwmrecordspt = prevwmrecords;
			wmrecordspt = wmrecords+1;
			wmrecords1pt = wmrecords;
			warpipt = warpi + 1;
			warpjpt = warpj + 1;
		}


		for( j=1; j<lastj; j++ )
		{
#ifdef xxxenablemultithread
//			fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
			if( chudanpt && *chudanpt != chudanref ) 
			{
//				fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
				*chudanres = 1;
				return( -1.0 );
			}
#endif
			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "\n i=%d, j=%d %c, %c", i, j, seq1[0][i], seq2[0][j] );
			fprintf( stderr, "%5.0f->", wm );
			fprintf( stderr, "%5.0f? (penal=%5.2f)", g=mi+*fgcp2pt*(1.0-gapfreq1[i]), *fgcp2pt*(1.0-gapfreq1[i]) );
#endif
			if( (g=mi+*fgcp2pt*gf1va) > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
//				fprintf( stderr, "Jump to %d (%c)!", mpi, seq2[0][mpi] );
			}
			if( (g=*prept+*ogcp2pt*gf1vapre) >= mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

#if 0 
			fprintf( stderr, "%5.0f->", wm );
			fprintf( stderr, "%5.0f? (penal=%5.2f)", g=*mjpt+fgcp1va*(1.0-gapfreq2[j]), fgcp1va*(1.0-gapfreq2[j]) );
#endif
			if( (g=*mjpt+ fgcp1va* *gf2pt) > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
//				fprintf( stderr, "Jump to %d (%c)!", *mpjpt, seq1[0][*mpjpt] );
			}
			if( (g=*prept+ ogcp1va* *gf2ptpre) >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif


			if( trywarp )
			{
#if USE_PENALTY_EX
				if( ( g=*prevwmrecordspt++ + fpenalty_shift + fpenalty_ex * ( i - prevwarpi[j-1] + j - prevwarpj[j-1] ) ) > wm ) // naka ha osokute kamawanai
#else
				if( ( g=*prevwmrecordspt++ + fpenalty_shift ) > wm ) // naka ha osokute kamawanai
#endif
				{
//					fprintf( stderr, "WARP!!\n" );
					if( warpn && prevwarpi[j-1] == warpis[warpn-1] && prevwarpj[j-1] == warpjs[warpn-1] )
					{
						*ijppt = warpbase + warpn - 1;
					}
					else
					{
						*ijppt = warpbase + warpn;
						warpis = realloc( warpis, sizeof(int) * ( warpn+1 ) );
						warpjs = realloc( warpjs, sizeof(int) * ( warpn+1 ) );
						warpis[warpn] = prevwarpi[j-1];
						warpjs[warpn] = prevwarpj[j-1];
						warpn++;
					}
					wm = g;
				}
	
#if 0
				fprintf( stderr, "%5.0f ", wm );
#endif
				curm = *curpt + wm;
	
				if( *wmrecords1pt > *wmrecordspt )
				{
					*wmrecordspt = *wmrecords1pt;
					*warpipt  = *(warpipt-1);
					*warpjpt  = *(warpjpt-1);
				}
				if( curm > *wmrecordspt )
				{
					*wmrecordspt = curm;
					*warpipt = i;
					*warpjpt = j;
				}
				wmrecordspt++;
				wmrecords1pt++;
				warpipt++;
				warpjpt++;
			}

			*curpt++ += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
			fgcp2pt++;
			ogcp2pt++;
			gf2ptpre++;
			gf2pt++;
		}
		lastverticalw[i] = currentw[lgth2-1];

		if( trywarp )
		{
			fltncpy( prevwmrecords, wmrecords, lastj );
			intncpy( prevwarpi, warpi, lastj );
			intncpy( prevwarpj, warpj, lastj );
		}
	}

	if( trywarp )
	{
//		fprintf( stderr, "wm = %f\n", wm );
//		fprintf( stderr, "warpn = %d\n", warpn );
		free( wmrecords );
		free( prevwmrecords );
		free( warpi );
		free( warpj );
		free( prevwarpi );
		free( prevwarpj );
	}

#if OUTGAP0TRY
	if( !outgap )
	{
		for( j=1; j<lgth2+1; j++ )
			currentw[j] -= offset * ( lgth2 - j ) / 2.0;
		for( i=1; i<lgth1+1; i++ )
			lastverticalw[i] -= offset * ( lgth1 - i  / 2.0);
	}
#endif
		
	/*
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr,"%s\n", seq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr,"%s\n", seq2[j] );
	fprintf( stderr, "====>" );
	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], seq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( mseq2[j], seq2[j] );
	*/
	if( localhom )
	{
		Atracking_localhom( impmatch, currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, warpis, warpjs, warpbase, &ngap1, &ngap2, reuseprofiles );
	}
	else
		Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, tailgp, warpis, warpjs, warpbase, &ngap1, &ngap2, reuseprofiles );

//	reporterr( "reuseprofiles after Atracking = %d\n", reuseprofiles );

	if( warpis ) free( warpis );
	if( warpjs ) free( warpjs );

//	fprintf( stderr, "### impmatch = %f\n", *impmatch );

	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}


	if( ngap1 || !reuseprofiles )
		for( i=0; i<icyc; i++ ) strcpy( seq1[i], mseq1[i] );
	if( ngap2 || !reuseprofiles )
		for( j=0; j<jcyc; j++ ) strcpy( seq2[j], mseq2[j] );

#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s\n", mseq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "%s\n", mseq2[j] );
#endif

//	fprintf( stdout, "firstmem=%d, icyc=%d, jcyc=%d, wm = %f\n", firstmem, icyc, jcyc, wm );

//	fprintf( stderr, "lgth1 = %d\n", lgth1 );
//	fprintf( stderr, "->      %d\n", strlen( seq1[0] ) );
	previousfirstlen = lgth1;
	previousfirstmem = firstmem;
	previousicyc = icyc;
	previouscall = calledbyfulltreebase;

	return( wm );
}

double A__align_gapmap( char **seq1, char **seq2, double *eff1, double *eff2, int icyc, int jcyc, int alloclen, LocalHom ***localhom, double *impmatch, int *gapmap1, int *gapmap2 )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{
	fprintf( stderr, "Unexpected error.  Please contact kazutaka.katoh@aist.go.jp\n" );
	exit( 1 );
}


double A__align_variousdist( int **which, double ***matrices, double **n_dynamicmtx, char **seq1, char **seq2, double *eff1, double *eff2, double **eff1s, double **eff2s, int icyc, int jcyc, int alloclen, LocalHom ***localhom, double *impmatch, char *sgap1, char *sgap2, char *egap1, char *egap2, int *chudanpt, int chudanref, int *chudanres, int headgp, int tailgp )
/* score no keisan no sai motokaraaru gap no atukai ni mondai ga aru */
{


//	int k;
	register int i, j, c;
	int ngap1, ngap2;
	int lasti, lastj;      /* outgap == 0 -> lgth1, outgap == 1 -> lgth1+1 */
	int lgth1, lgth2;
	int resultlen;
	double wm = 0.0;   /* int ?????? */
	double g;
	double *currentw, *previousw;
//	double fpenalty = (double)penalty;
#if USE_PENALTY_EX
	double fpenalty_ex = (double)penalty_ex;
#endif
#if 1
	double *wtmp;
	int *ijppt;
	double *mjpt, *prept, *curpt;
	int *mpjpt;
#endif
	static TLS double mi, *m;
	static TLS int **ijp;
	static TLS int mpi, *mp;
	static TLS double *w1, *w2;
	static TLS double *match;
	static TLS double *initverticalw;    /* kufuu sureba iranai */
	static TLS double *lastverticalw;    /* kufuu sureba iranai */
	static TLS char **mseq1;
	static TLS char **mseq2;
	static TLS char **mseq;
	static TLS double *ogcp1;
	static TLS double *ogcp2;
	static TLS double *fgcp1;
	static TLS double *fgcp2;
	static TLS double ***cpmx1s;
	static TLS double ***cpmx2s;
	static TLS int ***intwork;
	static TLS double ***doublework;
	static TLS int orlgth1 = 0, orlgth2 = 0;
	static TLS double *gapfreq1;
	static TLS double *gapfreq2;
	double fpenalty = (double)penalty;
	double fpenalty_shift = (double)penalty_shift;
	double *fgcp2pt;
	double *ogcp2pt;
	double fgcp1va;
	double ogcp1va;
	double *gf2pt;
	double *gf2ptpre;
	double gf1va;
	double gf1vapre;
	double headgapfreq1;
	double headgapfreq2;

	int *warpis = NULL;
	int *warpjs = NULL;
	int *warpi = NULL;
	int *warpj = NULL;
	int *prevwarpi = NULL;
	int *prevwarpj = NULL;
	double *wmrecords = NULL;
	double *prevwmrecords = NULL;
	int warpn = 0;
	int warpbase;
	double curm = 0.0;
	double *wmrecordspt, *wmrecords1pt, *prevwmrecordspt;
	int *warpipt, *warpjpt;
	int *nmask, **masklist1, **masklist2;


	if( seq1 == NULL )
	{
		if( orlgth1 )
		{
//			fprintf( stderr, "## Freeing local arrays in A__align\n" );
			orlgth1 = 0;
			orlgth2 = 0;

			imp_match_init_strict( NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, NULL );

			free( mseq1 );
			free( mseq2 );
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );

			FreeFloatVec( ogcp1 );
			FreeFloatVec( ogcp2 );
			FreeFloatVec( fgcp1 );
			FreeFloatVec( fgcp2 );


			FreeFloatCub( cpmx1s );
			FreeFloatCub( cpmx2s );

			FreeFloatVec( gapfreq1 );
			FreeFloatVec( gapfreq2 );

			FreeFloatCub( doublework );
			FreeIntCub( intwork );

		}
		else
		{
//			fprintf( stderr, "## Not allocated\n" );
		}
		return( 0.0 );
	}


#if SLOW
	nmask = calloc( maxdistclass, sizeof( int ) );
#else
	masklist1 = AllocateIntMtx( maxdistclass, 0 );
	masklist2 = AllocateIntMtx( maxdistclass, 0 );
	nmask = calloc( maxdistclass, sizeof( int ) );

	for( c=0; c<maxdistclass; c++ )
	{
		for( i=0; i<icyc; i++ ) for( j=0; j<jcyc; j++ )
		{
			if( eff1s[c][i] * eff2s[c][j] != 0.0 )
			{

				if( c != which[i][j] )
				{
					masklist1[c] = realloc( masklist1[c], sizeof( int ) * nmask[c]+1 );
					masklist2[c] = realloc( masklist2[c], sizeof( int ) * nmask[c]+1 );

					masklist1[c][nmask[c]] = i;
					masklist2[c][nmask[c]] = j;
					nmask[c]++;
				}
			}
		}
	}
	for( c=0; c<maxdistclass; c++ ) if( nmask[c] ) break;
	if( c<maxdistclass ) reporterr( "Found a complex grouping. This step may be a bit slow.\n" );
#endif

	lgth1 = strlen( seq1[0] );
	lgth2 = strlen( seq2[0] );
#if 0
	if( lgth1 == 0 || lgth2 == 0 )
	{
		fprintf( stderr, "WARNING (Aalignmm): lgth1=%d, lgth2=%d\n", lgth1, lgth2 );
	}
#endif
	if( lgth1 == 0 && lgth2 == 0 )
		return( 0.0 );

	if( lgth1 == 0 )
	{
		for( i=0; i<icyc; i++ )
		{
			j = lgth2;
			seq1[i][j] = 0;
			while( j ) seq1[i][--j] = *newgapstr;
//			fprintf( stderr, "seq1[i] = %s\n", seq1[i] );
		}
		return( 0.0 );
	}

	if( lgth2 == 0 )
	{
		for( i=0; i<jcyc; i++ )
		{
			j = lgth1;
			seq2[i][j] = 0;
			while( j ) seq2[i][--j] = *newgapstr;
//			fprintf( stderr, "seq2[i] = %s\n", seq2[i] );
		}
		return( 0.0 );
	}

	warpbase = lgth1 + lgth2;
	warpis = NULL;
	warpjs = NULL;
	warpn = 0;

	if( trywarp )
	{
//		fprintf( stderr, "In A__align_variousdist !!!!!\n" );
		if( headgp == 0 || tailgp == 0 )
		{
			fprintf( stderr, "At present, headgp and tailgp must be 1.\n" );
			exit( 1 );
		}
		wmrecords = AllocateFloatVec( lgth2+1 );
		warpi = AllocateIntVec( lgth2+1 );
		warpj = AllocateIntVec( lgth2+1 );
		prevwmrecords = AllocateFloatVec( lgth2+1 );
		prevwarpi = AllocateIntVec( lgth2+1 );
		prevwarpj = AllocateIntVec( lgth2+1 );
		for( i=0; i<lgth2+1; i++ ) wmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwmrecords[i] = 0.0;
		for( i=0; i<lgth2+1; i++ ) prevwarpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) prevwarpj[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpi[i] = -warpbase;
		for( i=0; i<lgth2+1; i++ ) warpj[i] = -warpbase;
	}

#if 0
	fprintf( stderr, "####  eff in SA+++align\n" );
	fprintf( stderr, "####  seq1[0] = %s\n", seq1[0] );
	fprintf( stderr, "####  strlen( seq1[0] ) = %d\n", strlen( seq1[0] ) );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "eff1[%d] = %f\n", i, eff1[i] );
	fprintf( stderr, "####  seq2[0] = %s\n", seq2[0] );
	fprintf( stderr, "####  strlen( seq2[0] ) = %d\n", strlen( seq2[0] ) );
	for( i=0; i<jcyc; i++ ) fprintf( stderr, "eff2[%d] = %f\n", i, eff2[i] );
#endif
	if( orlgth1 == 0 )
	{
		mseq1 = AllocateCharMtx( njob, 0 );
		mseq2 = AllocateCharMtx( njob, 0 );
	}



	if( lgth1 > orlgth1 || lgth2 > orlgth2 )
	{
		int ll1, ll2;


		if( orlgth1 > 0 && orlgth2 > 0 )
		{
			FreeFloatVec( w1 );
			FreeFloatVec( w2 );
			FreeFloatVec( match );
			FreeFloatVec( initverticalw );
			FreeFloatVec( lastverticalw );

			FreeFloatVec( m );
			FreeIntVec( mp );

			FreeCharMtx( mseq );

			FreeFloatVec( ogcp1 );
			FreeFloatVec( ogcp2 );
			FreeFloatVec( fgcp1 );
			FreeFloatVec( fgcp2 );


			FreeFloatCub( cpmx1s );
			FreeFloatCub( cpmx2s );

			FreeFloatVec( gapfreq1 );
			FreeFloatVec( gapfreq2 );

			FreeFloatCub( doublework );
			FreeIntCub( intwork );
		}

		ll1 = MAX( (int)(1.3*lgth1), orlgth1 ) + 100;
		ll2 = MAX( (int)(1.3*lgth2), orlgth2 ) + 100;

#if DEBUG
		fprintf( stderr, "\ntrying to allocate (%d+%d)xn matrices ... ", ll1, ll2 );
#endif

		w1 = AllocateFloatVec( ll2+2 );
		w2 = AllocateFloatVec( ll2+2 );
		match = AllocateFloatVec( ll2+2 );

		initverticalw = AllocateFloatVec( ll1+2 );
		lastverticalw = AllocateFloatVec( ll1+2 );

		m = AllocateFloatVec( ll2+2 );
		mp = AllocateIntVec( ll2+2 );

		mseq = AllocateCharMtx( njob, ll1+ll2 );

		ogcp1 = AllocateFloatVec( ll1+2 );
		ogcp2 = AllocateFloatVec( ll2+2 );
		fgcp1 = AllocateFloatVec( ll1+2 );
		fgcp2 = AllocateFloatVec( ll2+2 );

		cpmx1s = AllocateFloatCub( maxdistclass, nalphabets, ll1+2 );
		cpmx2s = AllocateFloatCub( maxdistclass, nalphabets, ll2+2 );

		gapfreq1 = AllocateFloatVec( ll1+2 );
		gapfreq2 = AllocateFloatVec( ll2+2 );

		doublework = AllocateFloatCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets ); 
		intwork = AllocateIntCub( maxdistclass, MAX( ll1, ll2 )+2, nalphabets+1 ); 

#if DEBUG
		fprintf( stderr, "succeeded\n" );
#endif

		orlgth1 = ll1 - 100;
		orlgth2 = ll2 - 100;
	}


	for( i=0; i<icyc; i++ )
	{
		mseq1[i] = mseq[i];
		seq1[i][lgth1] = 0;
	}
	for( j=0; j<jcyc; j++ )
	{
		mseq2[j] = mseq[icyc+j];
		seq2[j][lgth2] = 0;
	}


	if( orlgth1 > commonAlloc1 || orlgth2 > commonAlloc2 )
	{
		int ll1, ll2;

		if( commonAlloc1 && commonAlloc2 )
		{
			FreeIntMtx( commonIP );
		}

		ll1 = MAX( orlgth1, commonAlloc1 );
		ll2 = MAX( orlgth2, commonAlloc2 );

#if DEBUG
		fprintf( stderr, "\n\ntrying to allocate %dx%d matrices ... ", ll1+1, ll2+1 );
#endif

		commonIP = AllocateIntMtx( ll1+10, ll2+10 );

#if DEBUG
		fprintf( stderr, "succeeded\n\n" );
#endif

		commonAlloc1 = ll1;
		commonAlloc2 = ll2;
	}
	ijp = commonIP;

#if 0
	{
		double t = 0.0;
		for( i=0; i<icyc; i++ )
			t += eff1[i];
	fprintf( stderr, "## totaleff = %f\n", t );
	}
#endif

#if SLOW
#else
//	cpmx_calc_new( seq1, cpmx1, eff1, lgth1, icyc );
//	cpmx_calc_new( seq2, cpmx2, eff2, lgth2, jcyc );
	for( c=0; c<maxdistclass; c++ )
	{
		cpmx_calc_new( seq1, cpmx1s[c], eff1s[c], lgth1, icyc );
		cpmx_calc_new( seq2, cpmx2s[c], eff2s[c], lgth2, jcyc );
	}
#endif

	if( sgap1 )
	{
		new_OpeningGapCount( ogcp1, icyc, seq1, eff1, lgth1, sgap1 );
		new_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2, sgap2 );
		new_FinalGapCount( fgcp1, icyc, seq1, eff1, lgth1, egap1 );
		new_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2, egap2 );
		outgapcount( &headgapfreq1, icyc, sgap1, eff1 );
		outgapcount( &headgapfreq2, jcyc, sgap2, eff2 );
		outgapcount( gapfreq1+lgth1, icyc, egap1, eff1 );
		outgapcount( gapfreq2+lgth2, jcyc, egap2, eff2 );
	}
	else
	{
		st_OpeningGapCount( ogcp1, icyc, seq1, eff1, lgth1 );
		st_OpeningGapCount( ogcp2, jcyc, seq2, eff2, lgth2 );
		st_FinalGapCount( fgcp1, icyc, seq1, eff1, lgth1 );
		st_FinalGapCount( fgcp2, jcyc, seq2, eff2, lgth2 );
		headgapfreq1 = 0.0;
		headgapfreq2 = 0.0;
		gapfreq1[lgth1] = 0.0;
		gapfreq2[lgth2] = 0.0;
	}

	if( legacygapcost == 0 )
	{
		gapcountf( gapfreq1, seq1, icyc, eff1, lgth1 );
		gapcountf( gapfreq2, seq2, jcyc, eff2, lgth2 );
		for( i=0; i<lgth1+1; i++ ) gapfreq1[i] = 1.0 - gapfreq1[i];
		for( i=0; i<lgth2+1; i++ ) gapfreq2[i] = 1.0 - gapfreq2[i];
		headgapfreq1 = 1.0 - headgapfreq1;
		headgapfreq2 = 1.0 - headgapfreq2;
	}
	else
	{
		for( i=0; i<lgth1+1; i++ ) gapfreq1[i] = 1.0;
		for( i=0; i<lgth2+1; i++ ) gapfreq2[i] = 1.0;
		headgapfreq1 = 1.0;
		headgapfreq2 = 1.0;
	}

#if 0
	fprintf( stderr, "\ngapfreq1[] =" );
	for( i=0; i<lgth1; i++ ) fprintf( stderr, "%5.2f ", gapfreq1[i] );
	fprintf( stderr, "\n" );

	fprintf( stderr, "\ngapfreq2[] =" );
	for( i=0; i<lgth2; i++ ) fprintf( stderr, "%5.2f ", gapfreq2[i] );
	fprintf( stderr, "\n" );
#endif
	

	for( i=0; i<lgth1; i++ ) 
	{
		ogcp1[i] = 0.5 * ( 1.0 - ogcp1[i] ) * fpenalty * ( gapfreq1[i] );
		fgcp1[i] = 0.5 * ( 1.0 - fgcp1[i] ) * fpenalty * ( gapfreq1[i] );
	}

	for( i=0; i<lgth2; i++ ) 
	{
		ogcp2[i] = 0.5 * ( 1.0 - ogcp2[i] ) * fpenalty * ( gapfreq2[i] );
		fgcp2[i] = 0.5 * ( 1.0 - fgcp2[i] ) * fpenalty * ( gapfreq2[i] );
	}
#if 0
	for( i=0; i<lgth1; i++ ) 
		fprintf( stderr, "ogcp1[%d]=%f\n", i, ogcp1[i] );
#endif

	currentw = w1;
	previousw = w2;

//	for( i=0; i<icyc; i++ ) fprintf( stderr, "seq1[i] = %s\n", seq1[i] );
//	for( j=0; j<jcyc; j++ ) fprintf( stderr, "seq2[j] = %s\n", seq2[j] );

#if SLOW
	match_calc_slow( which, matrices, initverticalw, jcyc, seq2, eff2, icyc, seq1, eff1, 0, lgth1, *doublework, *intwork, 1, 1 );
#else
	fillzero( initverticalw, lgth1 );
	for( c=0; c<maxdistclass; c++ )
	{
//		fprintf( stderr, "c=%d matrices[c][W][W] = %f\n", c, matrices[c][amino_n['W']][amino_n['W']] );
//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "seq1[i] = %c, cpmx1s[c][3][%d] = %f\n", seq1[0][i], i, cpmx1s[c][3][i] );
//		for( i=0; i<lgth2; i++ ) fprintf( stderr, "seq2[i] = %c, cpmx2s[c][3][%d] = %f\n", seq2[0][i], i, cpmx2s[c][3][i] );
		match_calc_add( matrices[c], initverticalw, cpmx2s[c], cpmx1s[c], 0, lgth1, doublework[c], intwork[c], 1 );
//		for( i=0; i<lgth1; i++ ) fprintf( stderr, "c=%d, %d - %f\n", c, i, initverticalw[i] );
		if( nmask[c] ) match_calc_del( which, matrices, initverticalw, jcyc, seq2, eff2, icyc, seq1, eff1, 0, lgth1, c, nmask[c], masklist2[c], masklist1[c] );
	}
//	for( i=0; i<lgth1; i++ ) fprintf( stderr, "%d - %f\n", i, initverticalw[i] );
#endif

//	exit( 1 );

	if( localhom )
		imp_match_out_vead_tate( initverticalw, 0, lgth1 ); // 060306

#if SLOW
	match_calc_slow( which, matrices, currentw, icyc, seq1, eff1, jcyc, seq2, eff2, 0, lgth2, *doublework, *intwork, 1, 0 );
//	for( i=0; i<lgth2; i++ ) fprintf( stderr, "%d - %f\n", i, currentw[i] );
//	exit( 1 );
#else
	fillzero( currentw, lgth2 );
	for( c=0; c<maxdistclass; c++ )
	{
		match_calc_add( matrices[c], currentw, cpmx1s[c], cpmx2s[c], 0, lgth2, doublework[c], intwork[c], 1 );
		if( nmask[c] ) match_calc_del( which, matrices, currentw, icyc, seq1, eff1, jcyc, seq2, eff2, 0, lgth2, c, nmask[c], masklist1[c], masklist2[c] );
	}
//	for( i=0; i<lgth2; i++ ) fprintf( stderr, "%d - %f\n", i, currentw[i] );
//	exit( 1 );
#endif

	if( localhom )
		imp_match_out_vead( currentw, 0, lgth2 ); // 060306
#if 0 // -> tbfast.c
	if( localhom )
		imp_match_calc( n_dynamicmtx, currentw, icyc, jcyc, lgth1, lgth2, seq1, seq2, eff1, eff2, localhom, 1, 0 );

#endif

	if( headgp == 1 )
	{
		for( i=1; i<lgth1+1; i++ )
		{
//			initverticalw[i] += ( ogcp1[0] + fgcp1[i-1] ) ;
			initverticalw[i] += ( ogcp1[0] * headgapfreq2 + fgcp1[i-1] * gapfreq2[0] ) ;
		}
		for( j=1; j<lgth2+1; j++ )
		{
//			currentw[j] += ( ogcp2[0] + fgcp2[j-1] ) ;
			currentw[j] += ( ogcp2[0] * headgapfreq1 + fgcp2[j-1] * gapfreq1[0] ) ;
		}
	}
#if OUTGAP0TRY
	else
	{
		fprintf( stderr, "offset = %d\n", offset );
		for( j=1; j<lgth2+1; j++ )
			currentw[j] -= offset * j / 2.0;
		for( i=1; i<lgth1+1; i++ )
			initverticalw[i] -= offset * i / 2.0;
	}
#endif
#if 0
	fprintf( stderr, "\n   " );
	for( j=0; j<lgth2+1; j++ ) fprintf( stderr, "    %c ", seq2[0][j] );
	fprintf( stderr, "\n%c ", seq1[0][0]  );
	for( j=0; j<lgth2+1; j++ )
	{
		fprintf( stderr, "%5.0f ", currentw[j] );
	}
	fprintf( stderr, "\n"  );
#endif


	for( j=1; j<lgth2+1; ++j ) 
	{
//		m[j] = currentw[j-1] + ogcp1[1]; mp[j] = 0;
		m[j] = currentw[j-1] + ogcp1[1] * gapfreq2[j-1]; mp[j] = 0;;
	}
	if( lgth2 == 0 )
		lastverticalw[0] = 0.0; // Falign kara yobaretatoki kounarukanousei ari
	else
		lastverticalw[0] = currentw[lgth2-1];

	if( tailgp ) lasti = lgth1+1; else lasti = lgth1;
	lastj = lgth2+1;

#if XXXXXXX
fprintf( stderr, "currentw = \n" );
for( i=0; i<lgth1+1; i++ )
{
	fprintf( stderr, "%5.2f ", currentw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "initverticalw = \n" );
for( i=0; i<lgth2+1; i++ )
{
	fprintf( stderr, "%5.2f ", initverticalw[i] );
}
fprintf( stderr, "\n" );
fprintf( stderr, "fcgp\n" );
for( i=0; i<lgth1; i++ ) 
	fprintf( stderr, "fgcp1[%d]=%f\n", i, ogcp1[i] );
for( i=0; i<lgth2; i++ ) 
	fprintf( stderr, "fgcp2[%d]=%f\n", i, ogcp2[i] );
#endif

	for( i=1; i<lasti; i++ )
	{

#ifdef enablemultithread
//		fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
		if( chudanpt && *chudanpt != chudanref ) 
		{
//			fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
			if( masklist1 ) FreeIntMtx( masklist1 ); masklist1 = NULL;
			if( masklist2 ) FreeIntMtx( masklist2 ); masklist2 = NULL;
			if( nmask ) free( nmask ); nmask = NULL;
			*chudanres = 1;
			return( -1.0 );
		}
#endif
		wtmp = previousw; 
		previousw = currentw;
		currentw = wtmp;

		previousw[0] = initverticalw[i-1];

#if SLOW
		match_calc_slow( which, matrices, currentw, icyc, seq1, eff1, jcyc, seq2, eff2, i, lgth2, *doublework, *intwork, 0, 0 );
#else
		fillzero( currentw, lgth2 );
		for( c=0; c<maxdistclass; c++ )
		{
			match_calc_add( matrices[c], currentw, cpmx1s[c], cpmx2s[c], i, lgth2, doublework[c], intwork[c], 0 );
			if( nmask[c] ) match_calc_del( which, matrices, currentw, icyc, seq1, eff1, jcyc, seq2, eff2, i, lgth2, c, nmask[c], masklist1[c], masklist2[c] );
		}
#endif
#if 0
		if( i == 1 )
		{
			fprintf( stderr, "\n" );
			for( j=0; j<lgth2; j++ ) fprintf( stderr, "%d - %f\n", j, currentw[j] );
			exit( 1 );
		}
#endif

#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
		if( localhom )
		{
//			fprintf( stderr, "Calling imp_match_calc (o) lgth = %d, i = %d\n", lgth1, i );
#if  0
			imp_match_out_vead( currentw, i, lgth2 );
#else
			imp_match_out_vead( currentw, i, lgth2 );
#endif
		}
#if XXXXXXX
fprintf( stderr, "\n" );
fprintf( stderr, "i=%d\n", i );
fprintf( stderr, "currentw = \n" );
for( j=0; j<lgth2; j++ )
{
	fprintf( stderr, "%5.2f ", currentw[j] );
}
fprintf( stderr, "\n" );
#endif
		currentw[0] = initverticalw[i];

#if 0
		fprintf( stderr, "%c ", seq1[0][i] );
		for( j=0; j<lgth2+1; j++ )
		{
			fprintf( stderr, "%5.0f ", currentw[j] );
		}
		fprintf( stderr, "\n"  );
#endif
	
//		mi = previousw[0] + ogcp2[1]; mpi = 0;
		mi = previousw[0] + ogcp2[1] * gapfreq1[i-1]; mpi=0;
		ijppt = ijp[i] + 1;
		mjpt = m + 1;
		prept = previousw;
		curpt = currentw + 1;
		mpjpt = mp + 1;
		fgcp2pt = fgcp2;
		ogcp2pt = ogcp2 + 1;
		fgcp1va = fgcp1[i-1];
		ogcp1va = ogcp1[i];
		gf1va = gapfreq1[i];
		gf1vapre = gapfreq1[i-1];
		gf2pt = gapfreq2+1;
		gf2ptpre = gapfreq2;

		if( trywarp )
		{
			prevwmrecordspt = prevwmrecords;
			wmrecordspt = wmrecords+1;
			wmrecords1pt = wmrecords;
			warpipt = warpi + 1;
			warpjpt = warpj + 1;
		}

		for( j=1; j<lastj; j++ )
		{
#ifdef xxxenablemultithread
//			fprintf( stderr, "chudan = %d, %d\n", *chudanpt, chudanref );
			if( chudanpt && *chudanpt != chudanref ) 
			{
//				fprintf( stderr, "\n\n## CHUUDAN!!! S\n" );
				if( masklist1 ) FreeIntMtx( masklist1 ); masklist1 = NULL;
				if( masklist2 ) FreeIntMtx( masklist2 ); masklist2 = NULL;
				if( nmask ) free( nmask ); nmask = NULL;
				*chudanres = 1;
				return( -1.0 );
			}
#endif
			wm = *prept;
			*ijppt = 0;

#if 0
			fprintf( stderr, "\n i=%d, j=%d %c, %c", i, j, seq1[0][i], seq2[0][j] );
			fprintf( stderr, "%5.0f->", wm );
			fprintf( stderr, "%5.0f? (penal=%5.2f)", g=mi+*fgcp2pt*(1.0-gapfreq1[i]), *fgcp2pt*(1.0-gapfreq1[i]) );
#endif
			if( (g=mi+*fgcp2pt*gf1va) > wm )
			{
				wm = g;
				*ijppt = -( j - mpi );
//				fprintf( stderr, "Jump to %d (%c)!", mpi, seq2[0][mpi] );
			}
			if( (g=*prept+*ogcp2pt*gf1vapre) >= mi )
			{
				mi = g;
				mpi = j-1;
			}
#if USE_PENALTY_EX
			mi += fpenalty_ex;
#endif

#if 0 
			fprintf( stderr, "%5.0f->", wm );
			fprintf( stderr, "%5.0f? (penal=%5.2f)", g=*mjpt+fgcp1va*(1.0-gapfreq2[j]), fgcp1va*(1.0-gapfreq2[j]) );
#endif
			if( (g=*mjpt+ fgcp1va* *gf2pt) > wm )
			{
				wm = g;
				*ijppt = +( i - *mpjpt );
//				fprintf( stderr, "Jump to %d (%c)!", *mpjpt, seq1[0][*mpjpt] );
			}
			if( (g=*prept+ ogcp1va* *gf2ptpre) >= *mjpt )
			{
				*mjpt = g;
				*mpjpt = i-1;
			}
#if USE_PENALTY_EX
			m[j] += fpenalty_ex;
#endif
			if( trywarp )
			{
#if USE_PENALTY_EX
				if( ( g=*prevwmrecordspt++ + fpenalty_shift + fpenalty_ex * ( i - prevwarpi[j-1] + j - prevwarpj[j-1] ) ) > wm ) // naka ha osokute kamawanai
#else
				if( ( g=*prevwmrecordspt++ + fpenalty_shift ) > wm ) // naka ha osokute kamawanai
#endif
				{
					if( warpn && prevwarpi[j-1] == warpis[warpn-1] && prevwarpj[j-1] == warpjs[warpn-1] )
					{
						*ijppt = warpbase + warpn - 1;
					}
					else
					{
						*ijppt = warpbase + warpn;
						warpis = realloc( warpis, sizeof(int) * ( warpn+1 ) );
						warpjs = realloc( warpjs, sizeof(int) * ( warpn+1 ) );
						warpis[warpn] = prevwarpi[j-1];
						warpjs[warpn] = prevwarpj[j-1];
						warpn++;
					}
					wm = g;
				}
				curm = *curpt + wm;

				if( *wmrecords1pt > *wmrecordspt )
				{
					*wmrecordspt = *wmrecords1pt;
					*warpipt  = *(warpipt-1);
					*warpjpt  = *(warpjpt-1);
				}
				if( curm > *wmrecordspt )
				{
					*wmrecordspt = curm;
					*warpipt = i;
					*warpjpt = j;
				}
				wmrecordspt++;
				wmrecords1pt++;
				warpipt++;
				warpjpt++;
			}

#if 0
			fprintf( stderr, "%5.0f ", wm );
#endif
			*curpt++ += wm;
			ijppt++;
			mjpt++;
			prept++;
			mpjpt++;
			fgcp2pt++;
			ogcp2pt++;
			gf2ptpre++;
			gf2pt++;

		}
		lastverticalw[i] = currentw[lgth2-1];

		if( trywarp )
		{
			fltncpy( prevwmrecords, wmrecords, lastj );
			intncpy( prevwarpi, warpi, lastj );
			intncpy( prevwarpj, warpj, lastj );
		}
	}
	if( trywarp )
	{
//		fprintf( stderr, "wm = %f\n", wm );
//		fprintf( stderr, "warpn = %d\n", warpn );
		free( wmrecords );
		free( prevwmrecords );
		free( warpi );
		free( warpj );
		free( prevwarpi );
		free( prevwarpj );
	}


#if OUTGAP0TRY
	if( !outgap )
	{
		for( j=1; j<lgth2+1; j++ )
			currentw[j] -= offset * ( lgth2 - j ) / 2.0;
		for( i=1; i<lgth1+1; i++ )
			lastverticalw[i] -= offset * ( lgth1 - i  / 2.0);
	}
#endif
		
	/*
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr,"%s\n", seq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr,"%s\n", seq2[j] );
	fprintf( stderr, "====>" );
	for( i=0; i<icyc; i++ ) strcpy( mseq1[i], seq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( mseq2[j], seq2[j] );
	*/
	if( localhom )
	{
		Atracking_localhom( impmatch, currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, warpis, warpjs, warpbase, &ngap1, &ngap2, 0 );
	}
	else
		Atracking( currentw, lastverticalw, seq1, seq2, mseq1, mseq2, ijp, icyc, jcyc, tailgp, warpis, warpjs, warpbase, &ngap1, &ngap2, 0 );

	if( warpis ) free( warpis );
	if( warpjs ) free( warpjs );

//	fprintf( stderr, "### impmatch = %f\n", *impmatch );

	resultlen = strlen( mseq1[0] );
	if( alloclen < resultlen || resultlen > N )
	{
		fprintf( stderr, "alloclen=%d, resultlen=%d, N=%d\n", alloclen, resultlen, N );
		ErrorExit( "LENGTH OVER!\n" );
	}


	for( i=0; i<icyc; i++ ) strcpy( seq1[i], mseq1[i] );
	for( j=0; j<jcyc; j++ ) strcpy( seq2[j], mseq2[j] );
#if 0
	fprintf( stderr, "\n" );
	for( i=0; i<icyc; i++ ) fprintf( stderr, "%s\n", mseq1[i] );
	fprintf( stderr, "#####\n" );
	for( j=0; j<jcyc; j++ ) fprintf( stderr, "%s\n", mseq2[j] );
#endif

//	fprintf( stderr, "wm = %f\n", wm );


	if( masklist1 ) FreeIntMtx( masklist1 ); masklist1 = NULL;
	if( masklist2 ) FreeIntMtx( masklist2 ); masklist2 = NULL;
	if( nmask ) free( nmask ); nmask = NULL;

	return( wm );
}
