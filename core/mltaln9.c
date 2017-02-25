//tree methods and i think after group alignment tree construction
#include "mltaln.h"

#define DEBUG 0
#define CANONICALTREEFORMAT 1
#define MEMSAVE 1


//this method calculates the length of sequence without gaps
#if 0
int seqlen( char *seq )
{
	int val = 0;
	while( *seq )
		if( *seq++ != '-' ) val++;
	return( val );
}
#else
int seqlen( char *seq )
{
	int val = 0;
	if( *newgapstr == '-' )
	{
		while( *seq )
			if( *seq++ != '-' ) val++;
	}
	else
	{
		while( *seq )
		{
			if( *seq != '-' && *seq != *newgapstr ) val++;
			seq++;
		}
	}
	return( val );
}
#endif
//calculate length of an integer
int intlen( int *num )
{
	int value = 0;
	while( *num++ != -1 ) value++;
	return( value );
}
//check sequence characters and report error if unusual character is found
char seqcheck( char **seq )
{
	int i, len;
	char **seqbk = seq;
	while( *seq )	
	{
		len = strlen( *seq );
		for( i=0; i<len; i++ ) 
		{
			if( amino_n[(int)(*seq)[i]] == -1 ) 
			{

				reporterr(       "========================================================================= \n" );
				reporterr(       "========================================================================= \n" );
				reporterr(       "=== \n" );
				reporterr(       "=== Alphabet '%c' is unknown.\n", (*seq)[i] );
				reporterr(       "=== Please check site %d in sequence %d.\n", i+1, (int)(seq-seqbk+1) );
				reporterr(       "=== \n" );
				reporterr(       "=== To make an alignment having unusual characters (U, @, #, etc), try\n" );
				reporterr(       "=== %% mafft --anysymbol input > output\n" );
				reporterr(       "=== \n" );
				reporterr(       "========================================================================= \n" );
				reporterr(       "========================================================================= \n" );
				return( (int)(*seq)[i] );
			}
		}
		seq++;
	}
	return( 0 );
}
//concatenate two integers
void intcat( int *s1, int *s2 )
{
	while( *s1 != -1 ) s1++;
	while( *s2 != -1 ) 
	{
//		reporterr(       "copying %d\n", *s2 );
		*s1++ = *s2++;
	}
	*s1 = -1;
}
//copy integer s2 into s1
void intcpy( int *s1, int *s2 )
{
	while( *s2 != -1 ) 
	{
//		reporterr(       "copying %d\n", *s2 );
		*s1++ = *s2++;
	}
	*s1 = -1;
}
//copy n characters from integer s2 into s1
void intncpy( int *s1, int *s2, int n )
{
	while( n-- ) *s1++ = *s2++;
}
//copy n characters from float s2 into s1
void fltncpy( double *s1, double *s2, int n )
{
	while( n-- ) *s1++ = *s2++;
}
//count characters in an integer
static int countmem( int *s )
{
	int v = 0;
	while( *s++ != -1 ) v++;
	return( v );
}
//get last character in an integer
static int lastmem( int *s )
{
	while( *s++ != -1 ) 
		;
	return( *(s-2) );
}

	
void scmx_calc( int icyc, char **aseq, double *effarr, double **scmx )
{
	int  i, j, lgth;
	 
	lgth = strlen( aseq[0] );
	for( j=0; j<lgth; j++ )
	{
		for( i=0; i<nalphabets; i++ )
		{
			scmx[i][j] = 0;
		}
	}
	for( i=0; i<icyc+1; i++ )
	{
		int id;
		id = amino_n[(unsigned char)aseq[i][0]];
		scmx[id][0] += (double)effarr[i];
	}
	for( j=1; j<lgth-1; j++ )
	{
		for( i=0; i<icyc+1; i++ )
		{
			int id;
			id = amino_n[(unsigned char)aseq[i][j]];
			scmx[id][j] += (double)effarr[i];
		}
	}
	for( i=0; i<icyc+1; i++ )
	{
		int id;
		id = amino_n[(unsigned char)aseq[i][lgth-1]];
		scmx[id][lgth-1] += (double)effarr[i];
	}
}

void exitall( char arr[] )
{
	reporterr(       "%s\n", arr );
	exit( 1 );
}

void display( char **seq, int nseq )
{
	int i, imax;
	char b[121];

	if( !disp ) return;
		if( nseq > DISPSEQF ) imax = DISPSEQF;
		else                  imax = nseq;
		reporterr(       "    ....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+....,....+\n" );
		for( i=0; i<+imax; i++ )
		{
			strncpy( b, seq[i]+DISPSITEI, 120 );
			b[120] = 0;
			reporterr(       "%3d %s\n", i+1, b );
		}
}

void intergroup_score_consweight( char **seq1, char **seq2, double *eff1, double *eff2, int clus1, int clus2, int len, double *value )
{
	int i, j, k;
	int len2 = len - 2;
	unsigned char ms1, ms2;
	double tmpscore;
	char *mseq1, *mseq2;
	double efficient;

//	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
//	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];



	*value = 0.0;
	for( i=0; i<clus1; i++ ) 
	{
		for( j=0; j<clus2; j++ ) 
		{
			efficient = eff1[i] * eff2[j]; /* $B$J$<$+G[Ns$r;H$o$J$$$H$*$+$7$/$J$k(B, $BB?J,%P%0(B */
			mseq1 = seq1[i];
			mseq2 = seq2[j];
			tmpscore = 0.0;
			for( k=0; k<len; k++ ) 
			{
				ms1 = (unsigned char)mseq1[k];
				ms2 = (unsigned char)mseq2[k];
				if( ms1 == '-' && ms2 == '-' ) continue;
				tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
	
				if( ms1 == '-' ) 
				{
					tmpscore += (double)penalty;
					tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					while( (ms1=(unsigned char)mseq1[++k]) == '-' )
						;
//						tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					k--;
					if( k >len2 ) break;
					continue;
				}
				if( ms2 == '-' )
				{
					tmpscore += (double)penalty;
					tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					while( (ms2=(unsigned char)mseq2[++k]) == '-' )
						;
//						tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					k--;
					if( k > len2 ) break;
					continue;
				}
			}
			*value += (double)tmpscore * (double)efficient;
//			reporterr(       "val in _gapnomi = %f\n", *value );
		}
	}
#if 0
	fprintf( stdout, "###score = %f\n", score );
#endif
#if DEBUG
	reporterr(       "score in intergroup_score = %f\n", score );
#endif
//	return( score );
}
void intergroup_score_gapnomi( char **seq1, char **seq2, double *eff1, double *eff2, int clus1, int clus2, int len, double *value )
{
	int i, j, k;
	int len2 = len - 2;
	int ms1, ms2;
	double tmpscore;
	char *mseq1, *mseq2;
	double efficient;

//	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
//	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];



	*value = 0.0;
	for( i=0; i<clus1; i++ ) 
	{
		for( j=0; j<clus2; j++ ) 
		{
			efficient = eff1[i] * eff2[j]; /* $B$J$<$+G[Ns$r;H$o$J$$$H$*$+$7$/$J$k(B, $BB?J,%P%0(B */
			mseq1 = seq1[i];
			mseq2 = seq2[j];
			tmpscore = 0.0;
			for( k=0; k<len; k++ ) 
			{
				ms1 = (int)mseq1[k];
				ms2 = (int)mseq2[k];
				if( ms1 == (int)'-' && ms2 == (int)'-' ) continue;
//				tmpscore += (double)amino_dis[ms1][ms2];
	
				if( ms1 == (int)'-' ) 
				{
					tmpscore += (double)penalty;
//					tmpscore += (double)amino_dis[ms1][ms2];
					while( (ms1=(int)mseq1[++k]) == (int)'-' )
						;
//						tmpscore += (double)amino_dis[ms1][ms2];
					k--;
					if( k >len2 ) break;
					continue;
				}
				if( ms2 == (int)'-' )
				{
					tmpscore += (double)penalty;
//					tmpscore += (double)amino_dis[ms1][ms2];
					while( (ms2=(int)mseq2[++k]) == (int)'-' )
						;
//						tmpscore += (double)amino_dis[ms1][ms2];
					k--;
					if( k > len2 ) break;
					continue;
				}
			}
			*value += (double)tmpscore * (double)efficient;
//			reporterr(       "val in _gapnomi = %f\n", *value );
		}
	}
#if 0
	fprintf( stdout, "###score = %f\n", score );
#endif
#if DEBUG
	reporterr(       "score in intergroup_score = %f\n", score );
#endif
//	return( score );
}

void intergroup_score_multimtx( int **whichmtx, double ***scoringmatrices, char **seq1, char **seq2, double *eff1, double *eff2, int clus1, int clus2, int len, double *value )
{
	int i, j, k, c;
	int len2 = len - 2;
	int mn1, mn2;
	double tmpscore;
	char *mseq1, *mseq2;
	double efficient;
	int gapnum = amino_n['-'];

	double gaptmpscore;
	double gapscore = 0.0;

//	reporterr(       "#### in intergroup_score\n" );

//	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
//	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];

//	reporterr(       "\n intergroup_score_multimtx ..." );
	*value = 0.0;
	for( i=0; i<clus1; i++ ) 
	{
		for( j=0; j<clus2; j++ ) 
		{
			efficient = eff1[i] * eff2[j]; /* $B$J$<$+G[Ns$r;H$o$J$$$H$*$+$7$/$J$k(B, $BB?J,%P%0(B */
			c = whichmtx[i][j];
			mseq1 = seq1[i];
			mseq2 = seq2[j];
			tmpscore = 0.0;
			gaptmpscore = 0.0;
			for( k=0; k<len; k++ ) 
			{
				mn1 = amino_n[(unsigned char)(mseq1[k])];
				mn2 = amino_n[(unsigned char)(mseq2[k])];
				if( mn1 == gapnum && mn2 == gapnum ) continue;
				tmpscore += (double)scoringmatrices[c][mn1][mn2];
//				tmpscore += (double)scoringmtx[mn1][mn2];
	
				if( mn1 == gapnum ) 
				{
					tmpscore += (double)penalty;
					gaptmpscore += (double)penalty;
//					tmpscore += (double)scoringmtx[mn1][mn2];
					tmpscore += (double)scoringmatrices[c][mn1][mn2];
					while( (mn1=amino_n[(unsigned char)mseq1[++k]]) == gapnum )
						tmpscore += (double)scoringmatrices[c][mn1][mn2];
//						tmpscore += (double)scoringmtx[mn1][mn2];
					k--;
					if( k >len2 ) break;
					continue;
				}
				if( mn2 == gapnum )
				{
					tmpscore += (double)penalty;
					gaptmpscore += (double)penalty;
					tmpscore += (double)scoringmatrices[c][mn1][mn2];
//					tmpscore += (double)scoringmtx[mn1][mn2];
					while( (mn2=amino_n[(unsigned char)mseq2[++k]]) == gapnum )
						tmpscore += (double)scoringmatrices[c][mn1][mn2];
//						tmpscore += (double)scoringmtx[mn1][mn2];
					k--;
					if( k > len2 ) break;
					continue;
				}
			}
			*value += (double)tmpscore * (double)efficient;
			gapscore += (double)gaptmpscore * (double)efficient;
		}
	}
//	reporterr(       "done." );
#if 0
	reporterr(       "###gapscore = %f\n", gapscore );
#endif
#if DEBUG
	reporterr(       "score in intergroup_score = %f\n", score );
#endif
//	return( score );
}
void intergroup_score_dynmtx( double **offsetmtx, int scoringmtx[0x80][0x80], char **seq1, char **seq2, double *eff1, double *eff2, int clus1, int clus2, int len, double *value )
{
	int i, j, k;
	int len2 = len - 2;
	int ms1, ms2;
	double tmpscore;
	char *mseq1, *mseq2;
	double efficient;

	double gaptmpscore;
	double gapscore = 0.0;

//	reporterr(       "#### in intergroup_score\n" );

//	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
//	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];

	reporterr(       "\n intergroup_score_dynmtx ..." );
	*value = 0.0;
	for( i=0; i<clus1; i++ ) 
	{
		for( j=0; j<clus2; j++ ) 
		{
			efficient = eff1[i] * eff2[j]; /* $B$J$<$+G[Ns$r;H$o$J$$$H$*$+$7$/$J$k(B, $BB?J,%P%0(B */
			mseq1 = seq1[i];
			mseq2 = seq2[j];
			tmpscore = 0.0;
			gaptmpscore = 0.0;
			for( k=0; k<len; k++ ) 
			{
				ms1 = (int)mseq1[k];
				ms2 = (int)mseq2[k];
				if( ms1 == (int)'-' && ms2 == (int)'-' ) continue;
				tmpscore += (double)scoringmtx[ms1][ms2] + offsetmtx[i][j] * 600;
//				tmpscore += (double)scoringmtx[ms1][ms2];
	
				if( ms1 == (int)'-' ) 
				{
					tmpscore += (double)penalty;
					gaptmpscore += (double)penalty;
//					tmpscore += (double)scoringmtx[ms1][ms2];
					tmpscore += (double)scoringmtx[ms1][ms2] + offsetmtx[i][j] * 600;;
					while( (ms1=(int)mseq1[++k]) == (int)'-' )
						tmpscore += (double)scoringmtx[ms1][ms2] + offsetmtx[i][j] * 600;
//						tmpscore += (double)scoringmtx[ms1][ms2];
					k--;
					if( k >len2 ) break;
					continue;
				}
				if( ms2 == (int)'-' )
				{
					tmpscore += (double)penalty;
					gaptmpscore += (double)penalty;
					tmpscore += (double)scoringmtx[ms1][ms2] + offsetmtx[i][j] * 600;
//					tmpscore += (double)scoringmtx[ms1][ms2];
					while( (ms2=(int)mseq2[++k]) == (int)'-' )
						tmpscore += (double)scoringmtx[ms1][ms2] + offsetmtx[i][j] * 600;
//						tmpscore += (double)scoringmtx[ms1][ms2];
					k--;
					if( k > len2 ) break;
					continue;
				}
			}
			*value += (double)tmpscore * (double)efficient;
			gapscore += (double)gaptmpscore * (double)efficient;
		}
	}
	reporterr(       "done." );
#if 0
	reporterr(       "###gapscore = %f\n", gapscore );
#endif
#if DEBUG
	reporterr(       "score in intergroup_score = %f\n", score );
#endif
//	return( score );
}
void intergroup_score( char **seq1, char **seq2, double *eff1, double *eff2, int clus1, int clus2, int len, double *value )
{
	int i, j, k;
	int len2 = len - 2;
	unsigned char ms1, ms2;
	double tmpscore;
	char *mseq1, *mseq2;
	double efficient;

	double gaptmpscore;
	double gapscore = 0.0;

//	reporterr(       "#### in intergroup_score\n" );

//	totaleff1 = 0.0; for( i=0; i<clus1; i++ ) totaleff1 += eff1[i];
//	totaleff2 = 0.0; for( i=0; i<clus2; i++ ) totaleff2 += eff2[i];

	*value = 0.0;
	for( i=0; i<clus1; i++ ) 
	{
		for( j=0; j<clus2; j++ ) 
		{
			efficient = eff1[i] * eff2[j]; /* $B$J$<$+G[Ns$r;H$o$J$$$H$*$+$7$/$J$k(B, $BB?J,%P%0(B */
			mseq1 = seq1[i];
			mseq2 = seq2[j];
			tmpscore = 0.0;
			gaptmpscore = 0.0;
			for( k=0; k<len; k++ ) 
			{
				ms1 = (unsigned char)mseq1[k];
				ms2 = (unsigned char)mseq2[k];
				if( ms1 == '-' && ms2 == '-' ) continue;
//				tmpscore += (double)amino_dis[ms1][ms2];
				tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
	
				if( ms1 == '-' ) 
				{
					tmpscore += (double)penalty;
					gaptmpscore += (double)penalty;
//					tmpscore += (double)amino_dis[ms1][ms2];
					tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					while( (ms1=(unsigned char)mseq1[++k]) == '-' )
//						tmpscore += (double)amino_dis[ms1][ms2];
						tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					k--;
					if( k >len2 ) break;
					continue;
				}
				if( ms2 == '-' )
				{
					tmpscore += (double)penalty;
					gaptmpscore += (double)penalty;
//					tmpscore += (double)amino_dis[ms1][ms2];
					tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					while( (ms2=(unsigned char)mseq2[++k]) == '-' )
//						tmpscore += (double)amino_dis[ms1][ms2];
						tmpscore += (double)amino_dis_consweight_multi[ms1][ms2];
					k--;
					if( k > len2 ) break;
					continue;
				}
			}
			*value += (double)tmpscore * (double)efficient;
			gapscore += (double)gaptmpscore * (double)efficient;
		}
	}
#if 0
	reporterr(       "###gapscore = %f\n", gapscore );
#endif
#if DEBUG
	reporterr(       "score in intergroup_score = %f\n", score );
#endif
//	return( score );
}

double score_calc5( char **seq, int s, double **eff, int ex )  /* method 3 deha nai */
{
    int i, j, k;
    double c;
    int len = strlen( seq[0] );
    double score;
    double tmpscore;
    char *mseq1, *mseq2;
    double efficient;
#if DEBUG
	FILE *fp;
#endif

    score = 0.0;
    c = 0.0;

	for( i=0; i<s; i++ ) 
	{
		
			if( i == ex ) continue;
            efficient = eff[i][ex];
            mseq1 = seq[i];
            mseq2 = seq[ex];
            tmpscore = 0.0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(unsigned char)mseq1[k]][(unsigned char)mseq2[k]];

                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq1[++k] == '-' )
                        tmpscore += amino_dis[(unsigned char)mseq1[k]][(unsigned char)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq2[++k] == '-' )
                        tmpscore += amino_dis[(unsigned char)mseq1[k]][(unsigned char)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            score += (double)tmpscore * efficient;
/*
			fprintf( stdout, "%d-%d tmpscore = %f, eff = %f, tmpscore*eff = %f\n", i, ex, tmpscore, efficient, tmpscore*efficient );
*/
	}
	/*
	fprintf( stdout, "total score = %f\n", score );
	*/

    for( i=0; i<s-1; i++ )
    {
        for( j=i+1; j<s; j++ )
        {
			if( i == ex || j == ex ) continue;

            efficient = eff[i][j];
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0.0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(unsigned char)mseq1[k]][(unsigned char)mseq2[k]];

                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq1[++k] == '-' )
                        tmpscore += amino_dis[(unsigned char)mseq1[k]][(unsigned char)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq2[++k] == '-' )
                        tmpscore += amino_dis[(unsigned char)mseq1[k]][(unsigned char)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            score += (double)tmpscore * efficient;
        }
    }
/*
	reporterr(       "score in score_calc5 = %f\n", score );
*/
    return( (double)score );
/*

fprintf( trap_g, "score by fast = %f\n", (double)score );

tmpscore = score = 0.0;
	for( i=0; i<s; i++ ) 
	{
		if( i == ex ) continue;
		tmpscore = Cscore_m_1( seq, i, eff );
		fprintf( stdout, "%d %f\n", i, tmpscore );

		score += tmpscore;
	}
	tmpscore = Cscore_m_1( seq, ex, eff );
	fprintf( stdout, "ex%d %f\n", i, tmpscore );
	score += tmpscore;

	return( score );
*/
}


	
double score_calc4( char **seq, int s, double **eff, int ex )  /* method 3 deha nai */
{
    int i, j, k;
	double c;
    int len = strlen( seq[0] );
    double score;
    long tmpscore;
	char *mseq1, *mseq2;
	double efficient;

    score = 0.0;
	c = 0.0;
/*
	printf( "in score_calc4\n" );
	for( i=0; i<s; i++ )
	{
		for( j=0; j<s; j++ ) 
		{
			printf( "% 5.3f", eff[i][j] ); 
		}
		printf( "\n" );
		
	}
*/
    for( i=0; i<s-1; i++ )
    {
        for( j=i+1; j<s; j++ )
        {
			efficient = eff[i][j];
			if( mix == 1 ) efficient = 1.0;
			/*
			printf( "weight for %d v.s. %d = %f\n", i, j, efficient );
			*/
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(unsigned char)mseq1[k]][(unsigned char)mseq2[k]] + 400 * !scoremtx ;

				c += efficient;

                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty - n_dis[24][0];
                    while( mseq1[++k] == '-' )
						;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty - n_dis[24][0];
                    while( mseq2[++k] == '-' )
						;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
			/*
			if( x == 65 ) printf( "i=%d j=%d tmpscore=%d l=%d\n", i, j, tmpscore, len );
			*/
            score += (double)tmpscore * efficient;
        }
    }
    score /= c;
    return( (double)score );
}



void upg2( int nseq, double **eff, int ***topol, double **len )
{
    int i, j, k;
	double tmplen[M];

    static char **pair = NULL;

	if( !pair )
	{
		pair = AllocateCharMtx( njob, njob );
	}

	for( i=0; i<nseq; i++ ) tmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) pair[i][j] = 0;
    for( i=0; i<nseq; i++ ) pair[i][i] = 1;

    for( k=0; k<nseq-1; k++ )
    {
        double minscore = 9999.0;
        int im = -1, jm = -1;
        int count;

        for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
        {
            if( eff[i][j] < minscore )
            {
                minscore = eff[i][j];
                im = i; jm = j;
            }
        }
        for( i=0, count=0; i<nseq; i++ )
            if( pair[im][i] > 0 )
            {
                topol[k][0][count] = i;
                count++;
            }
        topol[k][0][count] = -1;
        for( i=0, count=0; i<nseq; i++ )
            if( pair[jm][i] > 0 )
            {
                topol[k][1][count] = i;
                count++;
            }
        topol[k][1][count] = -1;

		len[k][0] = minscore / 2.0 - tmplen[im];
		len[k][1] = minscore / 2.0 - tmplen[jm];

		tmplen[im] = minscore / 2.0;

        for( i=0; i<nseq; i++ ) pair[im][i] += ( pair[jm][i] > 0 );
        for( i=0; i<nseq; i++ ) pair[jm][i] = 0;

        for( i=0; i<nseq; i++ )
        {
            if( i != im && i != jm )
            {
                eff[MIN(i,im)][MAX(i,im)] =
                ( eff[MIN(i,im)][MAX(i,im)] + eff[MIN(i,jm)][MAX(i,jm)] ) / 2.0;
                eff[MIN(i,jm)][MAX(i,jm)] = 9999.0;
            }
            eff[im][jm] = 9999.0;
        }
#if DEBUG
        printf( "STEP-%03d:\n", k+1 );
		printf( "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) printf( " %03d", topol[k][0][i] );
        printf( "\n" );
		printf( "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) printf( " %03d", topol[k][1][i] );
        printf( "\n" );
#endif
    }
}

#define BLOCKSIZE 100
#define LARGEBLOCKSIZE 100

typedef struct _generaltdistarrthread_arg
{
	int para;
	int njob;
//	int thread_no;
	int m;
	int *nlen;
	char **seq;
	int **skiptable;
	int **pointt;
	int *ttable;
	int *tselfscore;
	int *posshared;
	int *joblist;
	double *result;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} generaldistarrthread_arg_t;

static void *generalkmerdistarrthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	generaldistarrthread_arg_t *targ = (generaldistarrthread_arg_t *)arg;
	int njob = targ->njob;
	int para = targ->para;
	int m = targ->m;
	int *nlen = targ->nlen;
	int **pointt = targ->pointt;
	int *ttable = targ->ttable;
	int *tselfscore = targ->tselfscore;
 	int *joblist = targ->joblist;
 	int *posshared = targ->posshared;
 	double *result = targ->result;
//	double **partmtx = targ->partmtx;
	int i, posinjoblist, n;

//			for( acpti=ac; acpti!=NULL; acpti=acpti->next )

		while( 1 )
		{
#ifdef enablemultithread
			if( para ) pthread_mutex_lock( targ->mutex );
#endif 
			if( *posshared >= njob ) // block no toki >=
			{
#ifdef enablemultithread
				if( para ) pthread_mutex_unlock( targ->mutex );
#endif
				commonsextet_p( NULL, NULL );
				return( NULL );
			}
			posinjoblist = *posshared;
			*posshared += LARGEBLOCKSIZE;
#ifdef enablemultithread
			if( para ) pthread_mutex_unlock( targ->mutex );
#endif

			for( n=0; n<LARGEBLOCKSIZE&&posinjoblist<njob; n++ )
			{
				i = joblist[posinjoblist++];
	
//				if( i == m ) continue; // iranai
		
				result[i] = distcompact( nlen[m], nlen[i], ttable, pointt[i], tselfscore[m], tselfscore[i] );
		
			}
		}
}

static void *generalmsadistarrthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	generaldistarrthread_arg_t *targ = (generaldistarrthread_arg_t *)arg;
	int njob = targ->njob;
	int para = targ->para;
	int m = targ->m;
	int *tselfscore = targ->tselfscore;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
 	int *joblist = targ->joblist;
 	int *posshared = targ->posshared;
 	double *result = targ->result;
//	double **partmtx = targ->partmtx;
	int i, posinjoblist, n;

//			for( acpti=ac; acpti!=NULL; acpti=acpti->next )

		while( 1 )
		{
#ifdef enablemultithread
			if( para ) pthread_mutex_lock( targ->mutex );
#endif 
			if( *posshared >= njob ) // block no toki >=
			{
#ifdef enablemultithread
				if( para ) pthread_mutex_unlock( targ->mutex );
#endif
				return( NULL );
			}
			posinjoblist = *posshared;
			*posshared += LARGEBLOCKSIZE;
#ifdef enablemultithread
			if( para ) pthread_mutex_unlock( targ->mutex );
#endif

			for( n=0; n<LARGEBLOCKSIZE&&posinjoblist<njob; n++ )
			{
				i = joblist[posinjoblist++];
	
//				if( i == m ) continue; // iranai
		
				result[i] = distcompact_msa( seq[m], seq[i], skiptable[m], skiptable[i], tselfscore[m], tselfscore[i] );
		
			}
		}
}

#if 1
static void kmerresetnearest( int nseq, Bchain *acpt, double **distfrompt, double *mindisfrompt, int *nearestpt, int pos, int *tselfscore, int **pointt, int *nlen, int *singlettable1, double *result, int *joblist )
{
	int i, j;
	double tmpdouble;
	double mindisfrom;
	int nearest;
//	double **effptpt;
	Bchain *acptj;
//	double *result;
//	int *joblist;

	mindisfrom = 999.9;
	nearest = -1;


//	reporterr( "resetnearest..\r" );
//	printf( "[%d], %f, dist=%d ->", pos, *distfrompt, *nearestpt );

//	mindisfrom = 999.9;
//	nearest = -1;


//	result = calloc( nseq, sizeof( double ) );
//	joblist = calloc( nseq, sizeof( int ) );


	for( acptj=(acpt+pos)->next,j=0; acptj!=NULL; acptj=acptj->next ) // setnearest ni awaseru
	{
		i = acptj->pos;
//		if( i == pos ) continue;

		if( distfrompt[pos] )
		{
			tmpdouble = result[i] = distfrompt[pos][i];
			if( tmpdouble < mindisfrom )
			{
				mindisfrom = tmpdouble;
				nearest = i;
			}
		}
		else if( distfrompt[i] )
		{
			tmpdouble = result[i] = distfrompt[i][pos];
			if( tmpdouble < mindisfrom )
			{
				mindisfrom = tmpdouble;
				nearest = i;
			}
		}
		else
			joblist[j++] = i;
	}

	for( acptj=acpt; (acptj&&acptj->pos!=pos); acptj=acptj->next ) // setnearest ni awaseru
	{
		i = acptj->pos;
//		if( i == pos ) continue;

		if( distfrompt[pos] )
		{
			tmpdouble = result[i] = distfrompt[pos][i];
			if( tmpdouble < mindisfrom )
			{
				mindisfrom = tmpdouble;
				nearest = i;
			}
		}
		else if( distfrompt[i] )
		{
			tmpdouble = result[i] = distfrompt[i][pos];
			if( tmpdouble < mindisfrom )
			{
				mindisfrom = tmpdouble;
				nearest = i;
			}
		}
		else
			joblist[j++] = i;
	}


	if( j )
	{
//		reporterr( "resetting in parallel!! j=%d\n", j );
//		exit( 1 );
		int posshared;
		generaldistarrthread_arg_t *targ;

#ifdef enablemultithread
		if( nthread )
		{
			pthread_t *handle;
			pthread_mutex_t mutex;

			targ = calloc( nthread, sizeof( generaldistarrthread_arg_t ) );
			handle = calloc( nthread, sizeof( pthread_t ) );
			posshared = 0;
			pthread_mutex_init( &mutex, NULL );
			for( i=0; i<nthread; i++ )
			{
				targ[i].para = 1;
				targ[i].njob = j;
				targ[i].m = pos;
				targ[i].tselfscore = tselfscore;
				targ[i].nlen = nlen;
				targ[i].pointt = pointt;
				targ[i].ttable = singlettable1;
				targ[i].joblist = joblist;
				targ[i].result = result;
				targ[i].posshared = &posshared;
				targ[i].mutex = &mutex;
	
				pthread_create( handle+i, NULL, generalkmerdistarrthread, (void *)(targ+i) );
			}
				
			for( j=0; j<nthread; j++ ) pthread_join( handle[j], NULL );
			pthread_mutex_destroy( &mutex );
			free( handle );
		}
		else
#endif
		{
			targ = calloc( 1, sizeof( generaldistarrthread_arg_t ) );
			posshared = 0;
			{
				targ[0].para = 0;
				targ[0].njob = j;
				targ[0].m = pos;
				targ[0].tselfscore = tselfscore;
				targ[0].nlen = nlen;
				targ[0].pointt = pointt;
				targ[0].ttable = singlettable1;
				targ[0].joblist = joblist;
				targ[0].result = result;
				targ[0].posshared = &posshared;
	
				generalkmerdistarrthread( targ );
			}
		}
		free( targ );
// sukoshi muda
		for( acptj=(acpt+pos)->next; acptj!=NULL; acptj=acptj->next ) // setnearest ni awaseru
		{
			j = acptj->pos;
			tmpdouble = result[j];
			if( tmpdouble < mindisfrom )
			{
				mindisfrom = tmpdouble;
				nearest = j;
			}
		}
	
		for( acptj=acpt; (acptj&&acptj->pos!=pos); acptj=acptj->next ) // setnearest ni awaseru
		{
			j = acptj->pos;
			tmpdouble = result[j];
			if( tmpdouble < mindisfrom )
			{
				mindisfrom = tmpdouble;
				nearest = j;
			}
		}
	}


	*mindisfrompt = mindisfrom;
	*nearestpt = nearest;

//	free( joblist );
//	free( result );
}

#else
static void kmerresetnearest( int nseq, Bchain *acpt, double **distfrompt, double *mindisfrompt, int *nearestpt, int pos, int *tselfscore, int **pointt, int *nlen, int *singlettable1, double *resultnotused, int *joblistnotused )
{
	int j;
	double tmpdouble;
	double mindisfrom;
	int nearest;
//	double **effptpt;
	Bchain *acptj;

	mindisfrom = 999.9;
	nearest = -1;


//	reporterr( "resetnearest..\r" );
//	printf( "[%d], %f, dist=%d ->", pos, *distfrompt, *nearestpt );

//	mindisfrom = 999.9;
//	nearest = -1;


	for( acptj=(acpt+pos)->next; acptj!=NULL; acptj=acptj->next ) // setnearest ni awaseru
	{
		j = acptj->pos;

		if( distfrompt[pos] )
			tmpdouble=distfrompt[pos][j];
		else if( distfrompt[j] )
			tmpdouble=distfrompt[j][pos];
//		else if( seq )
//			tmpdouble=distcompact_msa( seq[pos], seq[j], skiptable[pos], skiptable[j], tselfscore[pos], tselfscore[j] );
		else
			tmpdouble=distcompact( nlen[pos], nlen[j], singlettable1, pointt[j], tselfscore[pos], tselfscore[j] );
				

		if( tmpdouble < mindisfrom )
		{
			mindisfrom = tmpdouble;
			nearest = j;
		}
	}

	for( acptj=acpt; (acptj&&acptj->pos!=pos); acptj=acptj->next ) // setnearest ni awaseru
	{
		j = acptj->pos;

		if( distfrompt[pos] )
			tmpdouble=distfrompt[pos][j];
		else if( distfrompt[j] )
			tmpdouble=distfrompt[j][pos];
//		else if( seq )
//			tmpdouble=distcompact_msa( seq[pos], seq[j], skiptable[pos], skiptable[j], tselfscore[pos], tselfscore[j] );
		else
			tmpdouble=distcompact( nlen[pos], nlen[j], singlettable1, pointt[j], tselfscore[pos], tselfscore[j] );
				


		if( tmpdouble < mindisfrom )
		{
			mindisfrom = tmpdouble;
			nearest = j;
		}
	}
//	printf( "mindisfrom = %f\n", mindisfrom );

	*mindisfrompt = mindisfrom;
	*nearestpt = nearest;
}
#endif

#if 1
static void msaresetnearest( int nseq, Bchain *acpt, double **distfrompt, double *mindisfrompt, int *nearestpt, int pos, char **seq, int **skiptable, int *tselfscore, double *result, int *joblist )
{
	int i, j;
	double tmpdouble;
	double mindisfrom;
	int nearest;
//	double **effptpt;
	Bchain *acptj;
//	double *result;
//	int *joblist;

	mindisfrom = 999.9;
	nearest = -1;


//	reporterr( "resetnearest..\r" );
//	printf( "[%d], %f, dist=%d ->", pos, *distfrompt, *nearestpt );

//	mindisfrom = 999.9;
//	nearest = -1;


//	result = calloc( nseq, sizeof( double ) );
//	joblist = calloc( nseq, sizeof( int ) );

//	for( acptj=acpt,j=0; acptj!=NULL; acptj=acptj->next ) // setnearest ni awaseru
	for( acptj=(acpt+pos)->next,j=0; acptj!=NULL; acptj=acptj->next ) // setnearest ni awaseru
	{
		i = acptj->pos;
//		if( i == pos ) continue;

		if( distfrompt[pos] )
		{
			tmpdouble = result[i] = distfrompt[pos][i];
			if( tmpdouble < mindisfrom )
			{
				mindisfrom = tmpdouble;
				nearest = i;
			}
		}
		else if( distfrompt[i] )
		{
			tmpdouble = result[i] = distfrompt[i][pos];
			if( tmpdouble < mindisfrom )
			{
				mindisfrom = tmpdouble;
				nearest = i;
			}
		}
		else
			joblist[j++] = i;
	}

	for( acptj=acpt; (acptj&&acptj->pos!=pos); acptj=acptj->next ) // setnearest ni awaseru
	{
		i = acptj->pos;
//		if( i == pos ) continue;

		if( distfrompt[pos] )
		{
			tmpdouble = result[i] = distfrompt[pos][i];
			if( tmpdouble < mindisfrom )
			{
				mindisfrom = tmpdouble;
				nearest = i;
			}
		}
		else if( distfrompt[i] )
		{
			tmpdouble = result[i] = distfrompt[i][pos];
			if( tmpdouble < mindisfrom )
			{
				mindisfrom = tmpdouble;
				nearest = i;
			}
		}
		else
			joblist[j++] = i;
	}


	if( j )
	{
//		reporterr( "resetting in parallel!! j=%d\r", j );
//		exit( 1 );
		int posshared;
		generaldistarrthread_arg_t *targ;
		posshared = 0;

#ifdef enablemultithread
		if( nthread )
		{
			pthread_t *handle;
			pthread_mutex_t mutex;
			targ = calloc( nthread, sizeof( generaldistarrthread_arg_t ) );
			handle = calloc( nthread, sizeof( pthread_t ) );
			pthread_mutex_init( &mutex, NULL );
			for( i=0; i<nthread; i++ )
			{
				targ[i].para = 1;
				targ[i].njob = j;
				targ[i].m = pos;
				targ[i].tselfscore = tselfscore;
				targ[i].seq = seq;
				targ[i].skiptable = skiptable;
				targ[i].joblist = joblist;
				targ[i].result = result;
				targ[i].posshared = &posshared;
				targ[i].mutex = &mutex;
	
				pthread_create( handle+i, NULL, generalmsadistarrthread, (void *)(targ+i) );
			}
			for( j=0; j<nthread; j++ ) pthread_join( handle[j], NULL );
			pthread_mutex_destroy( &mutex );
			free( handle );
		}
		else
#endif
		{
			targ = calloc( 1, sizeof( generaldistarrthread_arg_t ) );
			{
				targ[0].para = 0;
				targ[0].njob = j;
				targ[0].m = pos;
				targ[0].tselfscore = tselfscore;
				targ[0].seq = seq;
				targ[0].skiptable = skiptable;
				targ[0].joblist = joblist;
				targ[0].result = result;
				targ[0].posshared = &posshared;
	
				generalmsadistarrthread( targ );
			}
		}
		free( targ );
// sukoshi muda
		for( acptj=(acpt+pos)->next; acptj!=NULL; acptj=acptj->next ) // setnearest ni awaseru
		{
			j = acptj->pos;
			tmpdouble = result[j];
			if( tmpdouble < mindisfrom )
			{
				mindisfrom = tmpdouble;
				nearest = j;
			}
		}
	
		for( acptj=acpt; (acptj&&acptj->pos!=pos); acptj=acptj->next ) // setnearest ni awaseru
		{
			j = acptj->pos;
			tmpdouble = result[j];
			if( tmpdouble < mindisfrom )
			{
				mindisfrom = tmpdouble;
				nearest = j;
			}
		}

	}

//	printf( "mindisfrom = %f\n", mindisfrom );

	*mindisfrompt = mindisfrom;
	*nearestpt = nearest;

//	free( joblist );
//	free( result );
}
#else
static void msaresetnearest( int nseq, Bchain *acpt, double **distfrompt, double *mindisfrompt, int *nearestpt, int pos, char **seq, int **skiptable, int *tselfscore, double *resultnotused, int *joblistnotused )
{
	int j;
	double tmpdouble;
	double mindisfrom;
	int nearest;
//	double **effptpt;
	Bchain *acptj;

	mindisfrom = 999.9;
	nearest = -1;


//	reporterr( "resetnearest..\r" );
//	printf( "[%d], %f, dist=%d ->", pos, *distfrompt, *nearestpt );

//	mindisfrom = 999.9;
//	nearest = -1;


	for( acptj=(acpt+pos)->next; acptj!=NULL; acptj=acptj->next ) // setnearest ni awaseru
	{
		j = acptj->pos;

		if( distfrompt[pos] )
			tmpdouble=distfrompt[pos][j];
		else if( distfrompt[j] )
			tmpdouble=distfrompt[j][pos];
		else
			tmpdouble=distcompact_msa( seq[pos], seq[j], skiptable[pos], skiptable[j], tselfscore[pos], tselfscore[j] );
//		else
//			tmpdouble=distcompact( nlen[pos], nlen[j], singlettable1, pointt[j], tselfscore[pos], tselfscore[j] );
				

		if( tmpdouble < mindisfrom )
		{
			mindisfrom = tmpdouble;
			nearest = j;
		}
	}

	for( acptj=acpt; (acptj&&acptj->pos!=pos); acptj=acptj->next ) // setnearest ni awaseru
	{
		j = acptj->pos;

		if( distfrompt[pos] )
			tmpdouble=distfrompt[pos][j];
		else if( distfrompt[j] )
			tmpdouble=distfrompt[j][pos];
		else
			tmpdouble=distcompact_msa( seq[pos], seq[j], skiptable[pos], skiptable[j], tselfscore[pos], tselfscore[j] );
//		else
//			tmpdouble=distcompact( nlen[pos], nlen[j], singlettable1, pointt[j], tselfscore[pos], tselfscore[j] );
				


		if( tmpdouble < mindisfrom )
		{
			mindisfrom = tmpdouble;
			nearest = j;
		}
	}
//	printf( "mindisfrom = %f\n", mindisfrom );

	*mindisfrompt = mindisfrom;
	*nearestpt = nearest;
}
#endif

static void setnearest( int nseq, Bchain *acpt, double **eff, double *mindisfrompt, int *nearestpt, int pos )
{
	int j;
	double tmpdouble;
	double mindisfrom;
	int nearest;
//	double **effptpt;
	Bchain *acptj;

	mindisfrom = 999.9;
	nearest = -1;

//	printf( "[%d], %f, dist=%d ->", pos, *mindisfrompt, *nearestpt );

//	if( (acpt+pos)->next ) effpt = eff[pos]+(acpt+pos)->next->pos-pos;

//	for( j=pos+1; j<nseq; j++ )
	for( acptj=(acpt+pos)->next; acptj!=NULL; acptj=acptj->next )
	{
		j = acptj->pos;
//		if( (tmpdouble=*effpt++) < *mindisfrompt )
		if( (tmpdouble=eff[pos][j-pos]) < mindisfrom )
		{
			mindisfrom = tmpdouble;
			nearest = j;
		}
	}
//	effptpt = eff;
//	for( j=0; j<pos; j++ )
	for( acptj=acpt; (acptj&&acptj->pos!=pos); acptj=acptj->next )
	{
		j = acptj->pos;
//		if( (tmpdouble=(*effptpt++)[pos-j]) < *mindisfrompt )
		if( (tmpdouble=eff[j][pos-j]) < mindisfrom )
		{
			mindisfrom = tmpdouble;
			nearest = j;
		}
	}

	*mindisfrompt = mindisfrom;
	*nearestpt = nearest;
//	printf( "%f, %d \n", pos, *mindisfrompt, *nearestpt );
}

static void setnearest_double_fullmtx( int nseq, Bchain *acpt, double **eff, double *mindisfrompt, int *nearestpt, int pos )
{
	int j;
	double tmpdouble;
	double **effptpt;
	Bchain *acptj;

	*mindisfrompt = 999.9;
	*nearestpt = -1;

//	if( (acpt+pos)->next ) effpt = eff[pos]+(acpt+pos)->next->pos-pos;

//	for( j=pos+1; j<nseq; j++ )
	for( acptj=(acpt+pos)->next; acptj!=NULL; acptj=acptj->next )
	{
		j = acptj->pos;
//		if( (tmpdouble=*effpt++) < *mindisfrompt )
		if( (tmpdouble=eff[pos][j]) < *mindisfrompt )
		{
			*mindisfrompt = tmpdouble;
			*nearestpt = j;
		}
	}
	effptpt = eff;
//	for( j=0; j<pos; j++ )
	for( acptj=acpt; (acptj&&acptj->pos!=pos); acptj=acptj->next )
	{
		j = acptj->pos;
//		if( (tmpdouble=(*effptpt++)[pos-j]) < *mindisfrompt )
		if( (tmpdouble=eff[j][pos]) < *mindisfrompt )
		{
			*mindisfrompt = tmpdouble;
			*nearestpt = j;
		}
	}
}



static void loadtreeoneline( int *ar, double *len, FILE *fp )
{
	static char gett[1000];
	int res;
	char *p;

	p = fgets( gett, 999, fp );
	if( p == NULL )
	{
		reporterr(       "\n\nFormat error (1) in the tree?  It has to be a bifurcated and rooted tree.\n" );
		reporterr(       "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n" );
		exit( 1 );
	}


	res = sscanf( gett, "%d %d %lf %lf", ar, ar+1, len, len+1 );
	if( res != 4 )
	{
		reporterr(       "\n\nFormat error (2) in the tree?  It has to be a bifurcated and rooted tree.\n" );
		reporterr(       "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n" );
		exit( 1 );
	}

	ar[0]--;
	ar[1]--;

	if( ar[0] >= ar[1] )
	{
		reporterr(       "\n\nIncorrect guide tree\n" );
		reporterr(       "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n" );
		exit( 1 );
	}


//	reporterr(       "ar[0] = %d, ar[1] = %d\n", ar[0], ar[1] );
//	reporterr(       "len[0] = %f, len[1] = %f\n", len[0], len[1] );
}

void loadtop( int nseq, double **mtx, int ***topol, double **len, char **name, int *nlen, Treedep *dep )
{
	int i, j, k, minijm, maxijm;
	int *intpt, *intpt2;
	int *hist = NULL;
	Bchain *ac = NULL;
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	int *pt1, *pt2, *pt11, *pt22;
	int *nmemar;
	int nmemim, nmemjm;
	char **tree;
	char *treetmp;
	char *nametmp, *nameptr, *tmpptr; 
	char namec;
	FILE *fp;
	int node[2];
	double *height;
	double clusterdist;
	int mpair, mi, mj;

	fp = fopen( "_guidetree", "r" );
	if( !fp )
	{
		reporterr(       "cannot open _guidetree\n" );
		exit( 1 );
	}

	if( !hist )
	{
		hist = AllocateIntVec( nseq );
		ac = (Bchain *)malloc( nseq * sizeof( Bchain ) );
		nmemar = AllocateIntVec( nseq );
//		treetmp = AllocateCharVec( nseq*50 );
		treetmp = NULL;
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( nseq, nseq*50 );
		tree = AllocateCharMtx( nseq, 0 );
		height = AllocateFloatVec( nseq );
	}

	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}


	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;


	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );
#if 0
		minscore = 999.9;
		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
		{
			i = acpti->pos;
//			reporterr(       "k=%d i=%d\n", k, i );
			if( mindisfrom[i] < minscore ) // muscle
			{
				im = i;
				minscore = mindisfrom[i];
			}
		}
		jm = nearest[im];
		if( jm < im ) 
		{
			j=jm; jm=im; im=j;
		}
#else
		len[k][0] = len[k][1] = -1.0;
		loadtreeoneline( node, len[k], fp );
		im = node[0];
		jm = node[1];

		if( im > nseq-1 || jm > nseq-1 || tree[im] == NULL || tree[jm] == NULL )
		{
			reporterr(       "\n\nCheck the guide tree.\n" );
			reporterr(       "im=%d, jm=%d\n", im+1, jm+1 );
			reporterr(       "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n" );
			exit( 1 );
		}

#endif

		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];

//		reporterr(       "prevnode = %d, nmemim = %d\n", prevnode, nmemim );

		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}


		nmemjm = nmemar[jm];
		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;

//		reporterr(       "prevnode = %d, nmemjm = %d\n", prevnode, nmemjm );

		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}


//		len[k][0] = ( minscore - tmptmplen[im] );
//		len[k][1] = ( minscore - tmptmplen[jm] );
//		len[k][0] = -1;
//		len[k][1] = -1;


		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;


		if( len[k][0] == -1 || len[k][1] == -1 )
		{
			reporterr( "Re-computing the length of branch %d..\n", k );
			clusterdist = 0.0;
			mpair = 0;
			for( i=0; (mi=topol[k][0][i])>-1; i++ ) for( j=0; (mj=topol[k][1][j])>-1; j++ ) 
			{
				minijm = MIN(mi,mj);
				maxijm = MAX(mi,mj);
				clusterdist += mtx[minijm][maxijm-minijm];
				mpair += 1;
			}
			clusterdist /= (double)mpair;
			reporterr( "clusterdist = %f\n", clusterdist );
			if( len[k][0] == -1 ) len[k][0] = clusterdist/2.0 - height[im];
			if( len[k][1] == -1 ) len[k][1] = clusterdist/2.0 - height[im];
	
			fprintf( stderr, "len0 = %f\n", len[k][0] );
			fprintf( stderr, "len1 = %f\n\n", len[k][1] );
		}

#if 0
        fprintf( stderr, "vSTEP-%03d:\n", k+1 );
		fprintf( stderr, "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) fprintf( stderr, " %03d", topol[k][0][i]+1 );
        fprintf( stderr, "\n" );
		fprintf( stderr, "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) fprintf( stderr, " %03d", topol[k][1][i]+1 );
        fprintf( stderr, "\n" );

#endif
		height[im] += len[k][0]; // for ig tree, 2015/Dec/25
		dep[k].distfromtip = height[im]; // for ig tree, 2015/Dec/25
//		reporterr( "##### dep[%d].distfromtip = %f\n", k, height[im] );



		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );

//		reporterr(       "im,jm=%d,%d\n", im, jm );
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
//		free( (void *)eff[jm] ); eff[jm] = NULL;

#if 0 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
//				reporterr(       "calling setnearest\n" );
//				setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#endif


    }
	fclose( fp );
	fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
		fprintf( fp, "#by loadtop\n" );
	fclose( fp );

	FreeCharMtx( tree );
	free( treetmp );
	free( nametmp );
	free( hist );
	free( (char *)ac );
	free( (void *)nmemar );
	free( height );

}


void stringshuffle( int *ary, int size )
{
	int i;
    for(i=0;i<size;i++)
    {
        int j = rand()%size;
        int t = ary[i];
        ary[i] = ary[j];
        ary[j] = t;
    }
}

void topolorder( int nseq, int *order, int *posinorder, int ***topol, Treedep *dep, int pos, int nchild )
{
	int *str;
	int child0, child1;
	child0 = dep[pos].child0;
	child1 = dep[pos].child1;

//	int i;


	if( nchild == 0 || nchild == 2 )
	{
		if( child0 == -1 )
		{
			str = calloc( 2, sizeof( int ) );
			str[0] = topol[pos][0][0]; // kanarazu memsave format nara, tanjunka dekiru.
			str[1] = -1;


//			for( i=0; order[i]!=-1; i++ )
//				;
//			reporterr( "0: i=%d, *posinorder=%d\n", i, *posinorder );

			intcpy( order+*posinorder, str );
//			intcat( order, str );


			*posinorder += 1;

			free( str );
		}
		else
		{
			topolorder( nseq, order, posinorder, topol, dep, child0, 2 );
		}
	}


	if( nchild == 1 || nchild == 2 )
	{
		if( child1 == -1 )
		{
			str = calloc( 2, sizeof( int ) );
			str[0] = topol[pos][1][0]; // kanarazu memsave format nara, tanjunka dekiru.
			str[1] = -1;


//			for( i=0; order[i]!=-1; i++ )
//				;
//			reporterr( "1: i=%d, *posinorder=%d\n", i, *posinorder );

			intcpy( order+*posinorder, str );
//			intcat( order, str );


			*posinorder += 1;
			free( str );
		}
		else
		{
			topolorder( nseq, order, posinorder, topol, dep, child1, 2 );
		}
	}
//	return( posinorder );
}

#if CANONICALTREEFORMAT
void createchain( int nseq, int ***topol, double **len, char **name, int *nlen, Treedep *dep, int treeout, int shuffle, int seed )
{
	FILE *fp;
	int i, j;
	double l, ll;
	int treelen;
	char **tree;
	char *instanttree;
	int posinit;
//	char *treetmp, *tt;
	char *nametmp, *nameptr, *tmpptr; 
	char namec;
	int *order;
	int im, jm, mm;

	if( treeout )
	{
//		treetmp = NULL;
		nametmp = AllocateCharVec( 1000 ); // nagasugi
		tree = AllocateCharMtx( nseq, 0 );

		treelen = nseq;
		for( i=0; i<nseq; i++ )
		{
	
			for( j=0; j<999; j++ ) nametmp[j] = 0;
			for( j=0; j<999; j++ ) 
			{
				namec = name[i][j];
				if( namec == 0 )
					break;
				else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
					nametmp[j] = namec;
				else
					nametmp[j] = '_';
			}
			nametmp[j] = 0;
//			sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
			if( outnumber )
				nameptr = strstr( nametmp, "_numo_e" ) + 8;
			else
				nameptr = nametmp + 1;
	
			if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame
	
			tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
			if( tree[i] == NULL )
			{
				reporterr(       "Cannot allocate tree!\n" );
				exit( 1 );
			}
			sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
			treelen += strlen( tree[i] ) + 20;

		}

		instanttree = calloc( treelen, sizeof( char ) );
		posinit = 0;
		for( i=0; i<nseq-1; i++ )
		{
			instanttree[i] = '(';
			posinit++;
		}

	}


	order = calloc( nseq, sizeof( int ) );
	for( i=0; i<nseq; i++ ) order[i] = i;

	srand( seed );
	if( shuffle ) stringshuffle( order, nseq );

	ll = l = 2.0 / nseq;

	im = order[0];
	jm = order[1];

	topol[0][0] = (int *)realloc( topol[0][0], ( 2 ) * sizeof( int ) );
	topol[0][1] = (int *)realloc( topol[0][1], ( 2 ) * sizeof( int ) );
	if( im < jm )
	{
		topol[0][0][0] = im;
		topol[0][0][1] = -1;
		topol[0][1][0] = jm;
		topol[0][1][1] = -1;
		mm = im;
	}	
	else
	{
		topol[0][0][0] = jm;
		topol[0][0][1] = -1;
		topol[0][1][0] = im;
		topol[0][1][1] = -1;
		mm = jm;
	}
	len[0][0] = len[0][1] = l;
	dep[0].child1 = -1;
	dep[0].child0 = -1;
	dep[0].distfromtip = l;
	ll += l;

	if( treeout )
	{
		posinit += sprintf( instanttree+posinit, "%s:%7.5f,", tree[im], len[0][0] );
//		reporterr( "instanttree = %s\n", instanttree );
	}

	for( i=1; i<nseq-1; i++ )
	{
		im = order[i];
		jm = order[i+1];

		if( mm < jm )
		{
#if MEMSAVE
			topol[i][0] = (int *)realloc( topol[i][0], ( 2 ) * sizeof( int ) );
			topol[i][0][0] = mm;
			topol[i][0][1] = -1;
#else
			topol[i][0] = (int *)realloc( topol[i][0], ( i + 2 ) * sizeof( int ) );
			intcpy( topol[i][0], topol[i-1][0] );
			intcat( topol[i][0], topol[i-1][1] );
#endif
			topol[i][1] = (int *)realloc( topol[i][1], ( 2 ) * sizeof( int ) );
			topol[i][1][0] = jm;
			topol[i][1][1] = -1;
	
//			reporterr( "step %d\n", i );
//			for( j=0; topol[i][0][j]!=-1; j++ ) reporterr( "%5d ", topol[i][0][j] );
//			reporterr( "\n", i );
//			for( j=0; topol[i][1][j]!=-1; j++ ) reporterr( "%5d ", topol[i][1][j] );
//			reporterr( "\n\n", i );
//			
			len[i][0] = l;
			len[i][1] = ll;
	
			if( dep ) 
			{
				dep[i].child0 = i-1;
				dep[i].child1 = -1;
				dep[i].distfromtip = ll;
			}
		}
		else
		{

#if MEMSAVE
			topol[i][1] = (int *)realloc( topol[i][1], ( 2 ) * sizeof( int ) );
			topol[i][1][0] = mm;
			topol[i][1][1] = -1;
#else
			topol[i][1] = (int *)realloc( topol[i][1], ( i + 2 ) * sizeof( int ) );
			intcpy( topol[i][1], topol[i-1][0] );
			intcat( topol[i][1], topol[i-1][1] );
#endif
			topol[i][0] = (int *)realloc( topol[i][0], ( 2 ) * sizeof( int ) );
			topol[i][0][0] = jm;
			topol[i][0][1] = -1;

			mm = jm;
	
//			reporterr( "step %d\n", i );
//			for( j=0; topol[i][0][j]!=-1; j++ ) reporterr( "%5d ", topol[i][0][j] );
//			reporterr( "\n", i );
//			for( j=0; topol[i][1][j]!=-1; j++ ) reporterr( "%5d ", topol[i][1][j] );
//			reporterr( "\n\n", i );
//			
	
			len[i][1] = l;
			len[i][0] = ll;
	
			if( dep ) 
			{
				dep[i].child1 = i-1;
				dep[i].child0 = -1;
				dep[i].distfromtip = ll;
			}
		}

		if( treeout ) 
		{
			posinit += sprintf( instanttree+posinit, "%s:%7.5f):%7.5f,", tree[im], ll-l, l );
//			reporterr( "instanttree (in loop) = %s\n", instanttree );
#if 0
			if( i % 1000 == 0 ) reporterr( "\r%d/%d", i, nseq );
//			reporterr( "size = %d\n", ( strlen( tree[im] ) + strlen( tree[jm] ) + 100 ) * sizeof( char ) );
//			reporterr( "size = %d\n", ( strlen( tree[im] ) + strlen( tree[jm] ) + 100 ) );
//			reporterr( "treetmp = %p\n", treetmp  );
			tt = realloc( treetmp, ( strlen( tree[im] ) + strlen( tree[jm] ) + 100 ) * sizeof( char ) ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
			if( tt == NULL )
			{
				reporterr(       "Cannot allocate treetmp\n" );
				exit( 1 );
			}
			treetmp = tt;
//			reporterr( "i=%d\n", i );
//			reporterr( "part1=%s\n", tree[0] );
//			reporterr( "part2=%s\n", tree[i+1] );
//			reporterr( "size = %d, %d\n", strlen( tree[0] ), strlen( tree[i+1] )  );
			sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[i][0], tree[jm], len[i][1] );
			free( tree[im] );
			free( tree[jm] );
			tree[jm] = calloc( strlen( treetmp )+1, sizeof( char ) );
			tree[im] = NULL;
			if( tree[jm] == NULL )
			{
				reporterr(       "Cannot reallocate tree!\n" );
				exit( 1 );
			}
			strcpy( tree[jm], treetmp );
#endif
		}
		ll += l;
	}
	if( treeout ) 
	{
		posinit += sprintf( instanttree+posinit, "%s:%7.5f)", tree[jm], ll-l );
		fp = fopen( "infile.tree", "w" );
//		fprintf( fp, "%s;\n", treetmp );
//		fprintf( fp, "#by createchain\n" );
		fprintf( fp, "%s;\n", instanttree );
		fclose( fp );
		FreeCharMtx( tree );
		free( nametmp );
		free( instanttree );
	}

	fp = fopen( "_guidetree", "w" );
	if( !fp )
	{
		reporterr(       "cannot open _guidetree\n" );
		exit( 1 );
	}
	for( i=0; i<nseq-1; i++ )
		fprintf( fp, "%d %d %f %f\n", topol[i][0][0]+1, topol[i][1][0]+1, len[i][0], len[i][1] );
	fclose( fp );

	free( order );

}
#else
void createchain( int nseq, int ***topol, double **len, char **name, int *nlen, Treedep *dep, int treeout, int shuffle, int seed )
{
	FILE *fp;
	int i, j;
	double l, ll;
	int treelen;
	char **tree;
	char *instanttree;
	int posinit;
//	char *treetmp, *tt;
	char *nametmp, *nameptr, *tmpptr; 
	char namec;
	int *order;
	int im, jm;

	if( treeout )
	{
//		treetmp = NULL;
		nametmp = AllocateCharVec( 1000 ); // nagasugi
		tree = AllocateCharMtx( nseq, 0 );

		treelen = nseq;
		for( i=0; i<nseq; i++ )
		{
	
			for( j=0; j<999; j++ ) nametmp[j] = 0;
			for( j=0; j<999; j++ ) 
			{
				namec = name[i][j];
				if( namec == 0 )
					break;
				else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
					nametmp[j] = namec;
				else
					nametmp[j] = '_';
			}
			nametmp[j] = 0;
//			sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
			if( outnumber )
				nameptr = strstr( nametmp, "_numo_e" ) + 8;
			else
				nameptr = nametmp + 1;
	
			if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame
	
			tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
			if( tree[i] == NULL )
			{
				reporterr(       "Cannot allocate tree!\n" );
				exit( 1 );
			}
			sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
			treelen += strlen( tree[i] ) + 20;

		}

		instanttree = calloc( treelen, sizeof( char ) );
		posinit = 0;
		for( i=0; i<nseq-1; i++ )
		{
			instanttree[i] = '(';
			posinit++;
		}

	}


	order = calloc( nseq, sizeof( int ) );
	for( i=0; i<nseq; i++ ) order[i] = i;

	srand( seed );
	if( shuffle ) stringshuffle( order, nseq );


	ll = l = 2.0 / nseq;

	for( i=0; i<nseq-1; i++ )
	{
		im = order[i];
		jm = order[i+1];

		topol[i][0] = (int *)realloc( topol[i][0], ( i + 2 ) * sizeof( int ) );
		topol[i][1] = (int *)realloc( topol[i][1], ( 2 ) * sizeof( int ) );

		for( j=0; j<i; j++ ) 
			topol[i][0][j] = order[j];
		topol[i][0][i] = im;
		topol[i][0][i+1] = -1;

		topol[i][1][0] = jm;
		topol[i][1][1] = -1;

		len[i][0] = l;
		len[i][1] = ll;
		ll += l;

		if( dep ) 
		{
			dep[i].child0 = i-1;
			dep[i].child1 = -1;
			dep[i].distfromtip = ll;
		}

		if( treeout ) 
		{
			if( i == 0 )
			{
				posinit += sprintf( instanttree+posinit, "%s:%7.5f,", tree[im], len[i][0] );
//				reporterr( "instanttree = %s\n", instanttree );
			}
			else if ( i == nseq-2 )
			{
				posinit += sprintf( instanttree+posinit, "%s:%7.5f):%7.5f,", tree[im], len[i-1][1], len[i-1][0] );
				posinit += sprintf( instanttree+posinit, "%s:%7.5f)", tree[jm], len[i][1] );
			}
			else
			{
				posinit += sprintf( instanttree+posinit, "%s:%7.5f):%7.5f,", tree[im], len[i-1][1], len[i-1][0] );
//				reporterr( "instanttree (in loop) = %s\n", instanttree );
			}
#if 0
			if( i % 1000 == 0 ) reporterr( "\r%d/%d", i, nseq );
//			reporterr( "size = %d\n", ( strlen( tree[im] ) + strlen( tree[jm] ) + 100 ) * sizeof( char ) );
//			reporterr( "size = %d\n", ( strlen( tree[im] ) + strlen( tree[jm] ) + 100 ) );
//			reporterr( "treetmp = %p\n", treetmp  );
			tt = realloc( treetmp, ( strlen( tree[im] ) + strlen( tree[jm] ) + 100 ) * sizeof( char ) ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
			if( tt == NULL )
			{
				reporterr(       "Cannot allocate treetmp\n" );
				exit( 1 );
			}
			treetmp = tt;
//			reporterr( "i=%d\n", i );
//			reporterr( "part1=%s\n", tree[0] );
//			reporterr( "part2=%s\n", tree[i+1] );
//			reporterr( "size = %d, %d\n", strlen( tree[0] ), strlen( tree[i+1] )  );
			sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[i][0], tree[jm], len[i][1] );
			free( tree[im] );
			free( tree[jm] );
			tree[jm] = calloc( strlen( treetmp )+1, sizeof( char ) );
			tree[im] = NULL;
			if( tree[jm] == NULL )
			{
				reporterr(       "Cannot reallocate tree!\n" );
				exit( 1 );
			}
			strcpy( tree[jm], treetmp );
#endif
		}
	}
	if( treeout ) 
	{
		fp = fopen( "infile.tree", "w" );
//		fprintf( fp, "%s;\n", treetmp );
//		fprintf( fp, "#by createchain\n" );
		fprintf( fp, "%s;\n", instanttree );
		fclose( fp );
		FreeCharMtx( tree );
		free( nametmp );
		free( instanttree );
	}

	fp = fopen( "_guidetree", "w" );
	if( !fp )
	{
		reporterr(       "cannot open _guidetree\n" );
		exit( 1 );
	}
#if CANONICALTREEFORMAT
	for( i=0; i<nseq-1; i++ )
		fprintf( fp, "%d %d %f %f\n", topol[i][0][0]+1, topol[i][1][0]+1, len[i][0], len[i][1] );
#else
	k = topol[0][0][0];
	for( i=0; i<nseq-1; i++ )
	{
		jm = topol[i][1][0];

		if( jm > k )
		{
			fprintf( fp, "%d %d %f %f\n", k+1, jm+1, len[i][0], len[i][1] );
		}
		else
		{
			fprintf( fp, "%d %d %f %f\n", jm+1, k+1, len[i][1], len[i][0] );
			k = jm;
		}
	}
#endif
	fclose( fp );
	free( order );
}
#endif

void loadtree( int nseq, int ***topol, double **len, char **name, int *nlen, Treedep *dep, int treeout )
{
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	int *hist = NULL;
	Bchain *ac = NULL;
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti;
	int *pt1, *pt2, *pt11, *pt22;
	int *nmemar;
	int nmemim, nmemjm;
	char **tree;
	char *treetmp;
	char *nametmp, *nameptr, *tmpptr; 
	char namec;
	FILE *fp;
	int node[2];
	double *height;

	fp = fopen( "_guidetree", "r" );
	if( !fp )
	{
		reporterr(       "cannot open _guidetree\n" );
		exit( 1 );
	}


	reporterr( "Loading a tree\n" );

	if( !hist )
	{
		hist = AllocateIntVec( nseq );
		ac = (Bchain *)malloc( nseq * sizeof( Bchain ) );
		nmemar = AllocateIntVec( nseq );
//		treetmp = AllocateCharVec( nseq*50 );
		if( dep ) height = AllocateFloatVec( nseq );
	}

	if( treeout )
	{
		treetmp = NULL;
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( nseq, nseq*50 );
		tree = AllocateCharMtx( nseq, 0 );

		for( i=0; i<nseq; i++ )
		{
			for( j=0; j<999; j++ ) nametmp[j] = 0;
			for( j=0; j<999; j++ ) 
			{
				namec = name[i][j];
				if( namec == 0 )
					break;
				else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
					nametmp[j] = namec;
				else
					nametmp[j] = '_';
			}
			nametmp[j] = 0;
//			sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
			if( outnumber )
				nameptr = strstr( nametmp, "_numo_e" ) + 8;
			else
				nameptr = nametmp + 1;
	
			if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame
	
			tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
			if( tree[i] == NULL )
			{
				reporterr(       "Cannot allocate tree!\n" );
				exit( 1 );
			}
			sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
		}

	}

	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;


	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );
#if 0
		minscore = 999.9;
		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
		{
			i = acpti->pos;
//			reporterr(       "k=%d i=%d\n", k, i );
			if( mindisfrom[i] < minscore ) // muscle
			{
				im = i;
				minscore = mindisfrom[i];
			}
		}
		jm = nearest[im];
		if( jm < im ) 
		{
			j=jm; jm=im; im=j;
		}
#else
		len[k][0] = len[k][1] = -1.0;
		loadtreeoneline( node, len[k], fp );
		im = node[0];
		jm = node[1];

//		if( im > nseq-1 || jm > nseq-1 || tree[im] == NULL || tree[jm] == NULL )
		if( im > nseq-1 || jm > nseq-1 )
		{
			reporterr(       "\n\nCheck the guide tree.\n" );
			reporterr(       "im=%d, jm=%d\n", im+1, jm+1 );
			reporterr(       "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n" );
			exit( 1 );
		}


		if( len[k][0] == -1.0 || len[k][1] == -1.0 )
		{
			reporterr(       "\n\nERROR: Branch length is not given.\n" );
			exit( 1 );
		}

		if( len[k][0] < 0.0 ) len[k][0] = 0.0;
		if( len[k][1] < 0.0 ) len[k][1] = 0.0;


#endif

		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];

//		reporterr(       "prevnode = %d, nmemim = %d\n", prevnode, nmemim );

		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}


		nmemjm = nmemar[jm];
		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;

//		reporterr(       "prevnode = %d, nmemjm = %d\n", prevnode, nmemjm );

		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}


//		len[k][0] = ( minscore - tmptmplen[im] );
//		len[k][1] = ( minscore - tmptmplen[jm] );
//		len[k][0] = -1;
//		len[k][1] = -1;


		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

//		mindisfrom[im] = 999.9;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
            }
        }


		if( treeout )
		{
			treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
			if( !treetmp )
			{
				reporterr(       "Cannot allocate treetmp\n" );
				exit( 1 );
			}
			sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
			free( tree[im] );
			free( tree[jm] );
			tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
			tree[jm] = NULL;
			if( tree[im] == NULL )
			{
				reporterr(       "Cannot reallocate tree!\n" );
				exit( 1 );
			}
			strcpy( tree[im], treetmp );
		}

//		reporterr(       "im,jm=%d,%d\n", im, jm );
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
//		free( (void *)eff[jm] ); eff[jm] = NULL;

#if 0 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
//				reporterr(       "calling setnearest\n" );
//				setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#endif


#if 0
        fprintf( stderr, "vSTEP-%03d:\n", k+1 );
		fprintf( stderr, "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) fprintf( stderr, " %03d", topol[k][0][i]+1 );
        fprintf( stderr, "\n" );
		fprintf( stderr, "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) fprintf( stderr, " %03d", topol[k][1][i]+1 );
        fprintf( stderr, "\n" );
#endif

		if( dep ) 
		{
			height[im] += len[k][0]; // for ig tree, 2015/Dec/25
			dep[k].distfromtip = height[im]; // for ig tree, 2015/Dec/25
//			reporterr(       "##### dep[%d].distfromtip = %f\n\n", k, height[im] );
		}

//		reporterr( "dep[%d].child0 = %d\n", k, dep[k].child0 );
//		reporterr( "dep[%d].child1 = %d\n", k, dep[k].child1 );
//		reporterr( "dep[%d].distfromtip = %f\n", k, dep[k].distfromtip );
    }
	fclose( fp );

	if( treeout )
	{
		fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s;\n", treetmp );
		fprintf( fp, "#by loadtree\n" );
		fclose( fp );
		FreeCharMtx( tree );
		free( treetmp );
		free( nametmp );
	}

	free( hist );
	free( (char *)ac );
	free( (void *)nmemar );
	if( dep ) free( height );

}

int check_guidetreefile( int *seed, int *npick, double *limitram )
{
	char string[100];
	char *sizestring;
	FILE *fp;
	double tanni;
	double tmpd;

	*seed = 0;
	*npick = 200;
	*limitram = 10.0 * 1000 * 1000 * 1000; // 10GB
	fp = fopen( "_guidetree", "r" );
	if( !fp )
	{
		reporterr(       "cannot open _guidetree\n" );
		exit( 1 );
	}

	fgets( string, 999, fp );
	fclose( fp );

	if( !strncmp( string, "shuffle", 7 ) )
	{
		sscanf( string+7, "%d", seed );
		reporterr( "shuffle, seed=%d\n", *seed );
		return( 's' );
	}
	else if( !strncmp( string, "pileup", 6 ) )
	{
		reporterr( "pileup.\n" );
		return( 'p' );
	}
	else if( !strncmp( string, "auto", 4 ) )
	{
		sscanf( string+4, "%d %d", seed, npick );
		reporterr( "auto, seed=%d, npick=%d\n", *seed, *npick );
		if( *npick < 2 )
		{
			reporterr( "Check npick\n" );
			exit( 1 );
		}
		return( 'a' );
	}
	else if( !strncmp( string, "test", 4 ) )
	{
		sscanf( string+4, "%d %d", seed, npick );
		reporterr( "calc, seed=%d, npick=%d\n", *seed, *npick );
		if( *npick < 2 )
		{
			reporterr( "Check npick\n" );
			exit( 1 );
		}
		return( 't' );
	}
	else if( !strncmp( string, "compact", 7 ) )
	{
		sizestring = string + 7;
		reporterr( "sizestring = %s\n", sizestring );
		if( strchr( sizestring, 'k' ) || strchr( sizestring, 'k' ) ) tanni = 1.0 * 1000; // kB
		else if( strchr( sizestring, 'M' ) || strchr( sizestring, 'm' ) ) tanni = 1.0 * 1000 * 1000; // GB
		else if( strchr( sizestring, 'G' ) || strchr( sizestring, 'g' ) ) tanni = 1.0 * 1000 * 1000 * 1000; // GB
		else if( strchr( sizestring, 'T' ) || strchr( sizestring, 't' ) ) tanni = 1.0 * 1000 * 1000 * 1000 * 1000; // TB
		else
		{
			reporterr( "\nSpecify initial ram usage by '--initialramusage xGB'\n\n\n" );
			exit( 1 );
		}
		sscanf( sizestring, "%lf", &tmpd );
		*limitram = tmpd * tanni;
		reporterr( "Initial RAM usage = %10.3fGB\n", *limitram/1000/1000/1000 );
		return( 'c' );
	}
	else if( !strncmp( string, "very compact", 12 ) )
	{
		reporterr( "very compact.\n" );
		return( 'C' );
	}
	else
	{
		reporterr( "loadtree.\n" );
		return( 'l' );
	}
}


static double sueff1, sueff05;
//static double sueff1_double, sueff05_double;

static double cluster_mix_double( double d1, double d2 )
{
	return( MIN( d1, d2 ) * sueff1 + ( d1 + d2 ) * sueff05 ); 
}
static double cluster_average_double( double d1, double d2 )
{
	return( ( d1 + d2 ) * 0.5 ); 
}
static double cluster_minimum_double( double d1, double d2 )
{
	return( MIN( d1, d2 ) ); 
}
#if 0
static double cluster_mix_double( double d1, double d2 )
{
	return( MIN( d1, d2 ) * sueff1_double + ( d1 + d2 ) * sueff05_double ); 
}
static double cluster_average_double( double d1, double d2 )
{
	return( ( d1 + d2 ) * 0.5 ); 
}
static double cluster_minimum_double( double d1, double d2 )
{
	return( MIN( d1, d2 ) ); 
}
#endif

static void increaseintergroupdistanceshalfmtx( double **eff, int ngroup, int **groups, int nseq )
{
	int nwarned = 0;
	int i, k, m, s1, s2, sl, ss;
	int *others, *tft;
	double maxdist, *dptr, dtmp;
	tft = calloc( nseq, sizeof( int * ) );
	others = calloc( nseq, sizeof( int * ) );

//	for( m=0; m<nseq-1; m++ ) for( k=m+1; k<nseq; k++ )
//		reporterr( "mtx[%d][%d] originally = %f (maxdist=%f)\n", m, k, eff[m][k-m], maxdist );

	reporterr( "\n" ); // Hitsuyou desu.
	for( i=0; i<ngroup; i++ )
	{
		if( groups[i][1] == -1 ) continue;

		for( m=0; m<nseq; m++ ) tft[m] = 0;
		for( m=0; (s1=groups[i][m])>-1; m++ ) tft[s1] = 1;
		for( m=0,k=0; m<nseq; m++ ) if( tft[m] == 0 ) others[k++] = m;
		others[k] = -1;

		maxdist = 0.0;
		for( m=1; (s2=groups[i][m])>-1; m++ ) for( k=0; (s1=groups[i][k])>-1&&k<m; k++ ) 
		{
//			reporterr( "m=%d, k=%d, s2=%d, s1=%d\n", m, k, s2, s1 );

			if( s2 > s1 )
			{
				sl = s2; ss = s1;
			}
			else
			{
				sl = s1; ss = s2;
			}
			dtmp = eff[ss][sl-ss];
			if( dtmp > maxdist ) maxdist = dtmp;
		}
//		reporterr( "maxdist = %f\n", maxdist );

		for( m=0; (s2=groups[i][m])>-1; m++ ) for( k=0; (s1=others[k])>-1; k++ ) 
		{
			if( s2 > s1 )
			{
				sl = s2; ss = s1;
			}
			else
			{
				sl = s1; ss = s2;
			}
			dptr = eff[ss] + sl-ss;
			if( *dptr < maxdist )
			{
				if( *dptr < 0.5 && nwarned++ < 100 ) reporterr( "# Sequences %d and %d seem to be closely related, but are not in the same sub MSA (%d) in your setting.\n", s2+1, s1+1, i+1 );
				*dptr = maxdist;
			}
		}
//		for( m=0; m<nseq-1; m++ ) for( k=m+1; k<nseq; k++ )
//			reporterr( "mtx[%d][%d] after modification%d = %f (maxdist=%f)\n", m, k, i, eff[m][k-m], maxdist );
	}
	if( nwarned > 100 ) reporterr( "# Sequenc.... (more pairs)\n" );

	free( tft );
	free( others );
}

static void increaseintergroupdistancesfullmtx( double **eff, int ngroup, int **groups, int nseq )
{
	int nwarned = 0;
	int i, k, m, s1, s2, sl, ss;
	int *others, *tft;
	double maxdist, *dptr, dtmp;
	tft = calloc( nseq, sizeof( int * ) );
	others = calloc( nseq, sizeof( int * ) );

	reporterr( "\n" ); // Hitsuyou desu.
	for( i=0; i<ngroup; i++ )
	{
		if( groups[i][1] == -1 ) continue;

		for( m=0; m<nseq; m++ ) tft[m] = 0;
		for( m=0; (s1=groups[i][m])>-1; m++ ) tft[s1] = 1;
		for( m=0,k=0; m<nseq; m++ ) if( tft[m] == 0 ) others[k++] = m;
		others[k] = -1;

		maxdist = 0.0;
		for( m=1; (s2=groups[i][m])>-1; m++ ) for( k=0; (s1=groups[i][k])>-1&&k<m; k++ ) 
		{
			if( s2 > s1 )
			{
				sl = s2; ss = s1;
			}
			else
			{
				sl = s1; ss = s2;
			}
			dtmp = eff[ss][sl];
			if( dtmp > maxdist ) maxdist = dtmp;
		}

//		reporterr( "maxdist = %f\n", maxdist );

		for( m=0; (s2=groups[i][m])>-1; m++ ) for( k=0; (s1=others[k])>-1; k++ ) 
		{
			if( s2 > s1 )
			{
				sl = s2; ss = s1;
			}
			else
			{
				sl = s1; ss = s2;
			}
			dptr = eff[ss] + sl;
			if( *dptr < maxdist )
			{
				if( *dptr < 0.5 && nwarned++ < 100 ) reporterr( "# Sequences %d and %d seem to be closely related, but are not in the same sub MSA (%d) in your setting.\n", s2+1, s1+1, i+1 );
				*dptr = maxdist;
			}
		}
	}
	if( nwarned > 100 ) reporterr( "# Sequenc.... (more pairs)\n" );

//	for( m=0; m<nseq-1; m++ ) for( k=m+1; k<nseq; k++ )
//		reporterr( "mtx[%d][%d] after modification = %f (maxdist=%f)\n", m, k, eff[m][k], maxdist );
	free( tft );
	free( others );
}

void fixed_supg_double_realloc_nobk_halfmtx_treeout_constrained( int nseq, double **eff, int ***topol, double **len, char **name, int *nlen, Treedep *dep, int ngroup, int **groups, int efffree )
{
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	double tmpdouble;
	double eff1, eff0;
	double *tmptmplen = NULL; //static?
	int *hist = NULL; //static?
	Bchain *ac = NULL; //static?
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti, *acptj;
	int *pt1, *pt2, *pt11, *pt22;
	int *nmemar; //static?
	int nmemim, nmemjm;
	double minscore;
	int *nearest = NULL; // by D.Mathog, a guess
	double *mindisfrom = NULL; // by D.Mathog, a guess
	char **tree; //static?
	char *treetmp; //static?
	char *nametmp, *nameptr, *tmpptr; //static?
	FILE *fp;
	double (*clusterfuncpt[1])(double,double);
	char namec;
	int *testtopol, **inconsistent;
	int **inconsistentpairlist;
	int ninconsistentpairs;
	int *warned;
	int allinconsistent;
	int firsttime;

	increaseintergroupdistanceshalfmtx( eff, ngroup, groups, nseq );

	sueff1 = 1 - (double)sueff_global;
	sueff05 = (double)sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_double;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_double;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_double;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}

	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateFloatVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateFloatVec( njob );
		nearest = AllocateIntVec( njob );
//		treetmp = AllocateCharVec( njob * ( B + 100 ) ); // nagasugi?
		treetmp = NULL; // kentou 2013/06/12
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( njob, njob*600 );
		tree = AllocateCharMtx( njob, 0 );
		testtopol = AllocateIntVec( njob + 1 );
		inconsistent = AllocateIntMtx( njob, njob ); // muda
		inconsistentpairlist = AllocateIntMtx( njob*(njob-1)/2+1, 2 ); // muda
		warned = AllocateIntVec( ngroup );
	}

	
	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

	for( i=0; i<nseq; i++ ) setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	ninconsistentpairs = 0;
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );

		for( i=0; i<ninconsistentpairs; i++ ) inconsistent[inconsistentpairlist[i][0]][inconsistentpairlist[i][1]] = 0;
//		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) for( acptj=acpti->next; acptj!=NULL; acptj=acptj->next ) inconsistent[acpti->pos][acptj->pos] = 0; // osoi!!!
		ninconsistentpairs = 0;
		firsttime = 1;
		while( 1 )
		{
			if( firsttime )
			{
				firsttime = 0;
				minscore = 999.9;
				for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
				{
					i = acpti->pos;
//					reporterr(       "k=%d i=%d\n", k, i );
					if( mindisfrom[i] < minscore ) // muscle
					{
						im = i;
						minscore = mindisfrom[i];
					}
				}
				jm = nearest[im];
				if( jm < im ) 
				{
					j=jm; jm=im; im=j;
				}
			}
			else
			{
				minscore = 999.9;
				for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
				{
					i = acpti->pos;
//					reporterr(       "k=%d i=%d\n", k, i );
					for( acptj=acpti->next; acptj!=NULL; acptj=acptj->next )
					{
						j = acptj->pos;
						if( !inconsistent[i][j] && (tmpdouble=eff[i][j-i]) < minscore )
						{
							minscore = tmpdouble;
							im = i; jm = j;
						}
					}
					for( acptj=ac; (acptj&&acptj->pos!=i); acptj=acptj->next )
					{
						j = acptj->pos;
						if( !inconsistent[j][i] && (tmpdouble=eff[j][i-j]) < minscore )
						{
							minscore = tmpdouble;
							im = j; jm = i;
						}
					}
				}
			}


			allinconsistent = 1;
			for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
			{
				for( acptj=acpti->next; acptj!=NULL; acptj=acptj->next )
				{
					if( inconsistent[acpti->pos][acptj->pos] == 0 )
					{
						allinconsistent = 0;
						goto exitloop_f;
					}
				}
			}
			exitloop_f:

			if( allinconsistent )
			{
				reporterr(       "\n\n\nPlease check whether the grouping is possible.\n\n\n" );
				exit( 1 );
			}
#if 1
			intpt = testtopol;
			prevnode = hist[im];
			if( prevnode == -1 )
			{
				*intpt++ = im;
			}
			else
			{
				for( intpt2=topol[prevnode][0]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
				for( intpt2=topol[prevnode][1]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
			}
			
			prevnode = hist[jm];
			if( prevnode == -1 )
			{
				*intpt++ = jm;
			}
			else
			{
				for( intpt2=topol[prevnode][0]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
				for( intpt2=topol[prevnode][1]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
			}
			*intpt = -1;
//			reporterr(       "testtopol = \n" );
//       	for( i=0; testtopol[i]>-1; i++ ) reporterr(       " %03d", testtopol[i]+1 );
//			reporterr(       "\n" );
#endif
			for( i=0; i<ngroup; i++ )
			{
//				reporterr(       "groups[%d] = \n", i );
//  		     	for( j=0; groups[i][j]>-1; j++ ) reporterr(       " %03d", groups[i][j]+1 );
//				reporterr(       "\n" );
				if( overlapmember( groups[i], testtopol ) )
				{
					if( !includemember( testtopol, groups[i] ) && !includemember( groups[i], testtopol ) )
					{
						if( !warned[i] )
						{
							warned[i] = 1;
							reporterr(       "\n###################################################################\n" );
							reporterr(       "# WARNING: Group %d is forced to be a monophyletic cluster.\n", i+1 );
							reporterr(       "###################################################################\n" );
						}
						inconsistent[im][jm] = 1;
						inconsistentpairlist[ninconsistentpairs][0] = im;
						inconsistentpairlist[ninconsistentpairs][1] = jm;
						ninconsistentpairs++;
						break;
					}
				}
			}
			if( i == ngroup )
			{
//				reporterr(       "OK\n" );
				break;
			}
		}


		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];
		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;
		nmemjm = nmemar[jm];
		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );
		if( len[k][0] < 0.0 ) len[k][0] = 0.0;
		if( len[k][1] < 0.0 ) len[k][1] = 0.0;

		if( dep ) dep[k].distfromtip = minscore;
//		reporterr(       "\n##### dep[%d].distfromtip = %f\n", k, minscore );

		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

		mindisfrom[im] = 999.9;
		eff[im][jm-im] = 999.9;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim-miniim];
				eff1 = eff[minijm][maxijm-minijm];
#if 0
                		tmpdouble = eff[miniim][maxiim-miniim] =
				MIN( eff0, eff1 ) * sueff1 + ( eff0 + eff1 ) * sueff05; 
#else
                tmpdouble = eff[miniim][maxiim-miniim] =
				(clusterfuncpt[0])( eff0, eff1 );
#endif
#if 1
				if( tmpdouble < mindisfrom[i]  )
				{
					mindisfrom[i] = tmpdouble;
					nearest[i] = im;
				}
				if( tmpdouble < mindisfrom[im]  )
				{
					mindisfrom[im] = tmpdouble;
					nearest[im] = i;
				}
				if( nearest[i] == jm )
				{
					nearest[i] = im;
				}
#endif
            }
        }

		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );

		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
		if( efffree )
		{
			free( (void *)eff[jm] ); eff[jm] = NULL;
		}

#if 1 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
				}
				if( eff[miniim][maxiim-miniim] > mindisfrom[i] )
					setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#endif


#if 0
        reporterr(       "\noSTEP-%03d:\n", k+1 );
		reporterr(       "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) reporterr(       " %03d", topol[k][0][i]+1 );
        reporterr(       "\n" );
		reporterr(       "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) reporterr(       " %03d", topol[k][1][i]+1 );
        reporterr(       "\n\n" );
#endif
    }
	fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
	fclose( fp );

	free( tree[0] );
	free( tree );
	free( treetmp );
	free( nametmp );
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	free( mindisfrom );
	free( nearest );
	free( testtopol );
	FreeIntMtx( inconsistent );
	FreeIntMtx( inconsistentpairlist );
	free( warned );
}

void makecompositiontable_global( int *table, int *pointt )
{
	int point;

	while( ( point = *pointt++ ) != END_OF_VEC )
		table[point]++;
}

typedef struct _resetnearestthread_arg
{
	int para;
//	int thread_no;
	int im;
	int nseq;
	double **partmtx;
	double *mindist;
	int *nearest;
	char **seq;
	int **skiptable;
	int *tselfscore;
	int **pointt;
	int *nlen;
	double *result;
	int *joblist;
 	Bchain **acpt;
 	Bchain *ac;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} resetnearestthread_arg_t;

static void *msaresetnearestthread( void *arg )
{
	resetnearestthread_arg_t *targ = (resetnearestthread_arg_t *)arg;
//	int thread_no = targ->thread_no;
	int para = targ->para;
	int im = targ->im;
	int nseq = targ->nseq;
	double **partmtx = targ->partmtx;
	double *mindist = targ->mindist;
	int *nearest = targ->nearest;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	int *tselfscore = targ->tselfscore;
	double *result = targ->result;
	int *joblist = targ->joblist;
 	Bchain **acpt = targ->acpt;
 	Bchain *ac = targ->ac;

	Bchain *acptbk;
	Bchain *acptinit;

	int i;

	acptinit = *acpt;
	while( 1 )
	{
#ifdef enablemultithread
		if( para ) pthread_mutex_lock( targ->mutex );
#endif 
		if( *acpt == NULL )
		{
#ifdef enablemultithread
			if( para ) pthread_mutex_unlock( targ->mutex );
#endif
			commonsextet_p( NULL, NULL );
			return( NULL );
		}
		acptbk = *acpt;
		*acpt = (*acpt)->next;

#ifdef enablemultithread
		if( para ) pthread_mutex_unlock( targ->mutex );
#endif
		i = acptbk->pos;
		if( nearest[i] == im ) 
		{
			if( partmtx[im][i] > mindist[i] )
			{
				msaresetnearest( nseq, ac, partmtx, mindist+i, nearest+i, i, seq, skiptable, tselfscore, result, joblist );
			}
		}
	}
}

static void *kmerresetnearestthread( void *arg )
{
	resetnearestthread_arg_t *targ = (resetnearestthread_arg_t *)arg;
//	int thread_no = targ->thread_no;
	int para = targ->para;
	int im = targ->im;
	int nseq = targ->nseq;
	double **partmtx = targ->partmtx;
	double *mindist = targ->mindist;
	int *nearest = targ->nearest;
	int *tselfscore = targ->tselfscore;
	int **pointt = targ->pointt;
	int *nlen = targ->nlen;
	double *result = targ->result;
	int *joblist = targ->joblist;
 	Bchain **acpt = targ->acpt;
 	Bchain *ac = targ->ac;

	int *singlettable1;

	Bchain *acptbk;
	Bchain *acptinit;

	int i;

	acptinit = *acpt;
	while( 1 )
	{
#ifdef enablemultithread
		if( para ) pthread_mutex_lock( targ->mutex );
#endif 
		if( *acpt == NULL )
		{
#ifdef enablemultithread
			if( para ) pthread_mutex_unlock( targ->mutex );
#endif
			commonsextet_p( NULL, NULL );
			return( NULL );
		}
		acptbk = *acpt;
		*acpt = (*acpt)->next;

#ifdef enablemultithread
		if( para ) pthread_mutex_unlock( targ->mutex );
#endif
		i = acptbk->pos;
		if( nearest[i] == im ) 
		{
			if( partmtx[im][i] > mindist[i] )
			{
				if( pointt ) // kmer
				{
					singlettable1 = (int *)calloc( tsize, sizeof( int ) );
					makecompositiontable_global( singlettable1, pointt[i] );
				}
				kmerresetnearest( nseq, ac, partmtx, mindist+i, nearest+i, i, tselfscore, pointt, nlen, singlettable1, result, joblist );
				if( pointt ) free( singlettable1 ); singlettable1 = NULL;// kmer
				if( pointt ) commonsextet_p( NULL, NULL );
			}
		}
	}
}


typedef struct _compactdistarrthread_arg
{
	int para;
	int njob;
//	int thread_no;
	int im;
	int jm;
	int *nlen;
	char **seq;
	int **skiptable;
	int **pointt;
	int *table1;
	int *table2;
	int *tselfscore;
	Bchain **acpt;
	int *posshared;
	double *mindist;
	double *newarr;
	double **partmtx;
	int *nearest;
	int *joblist;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} compactdistarrthread_arg_t;

static void *verycompactkmerdistarrthreadjoblist( void *arg ) // enablemultithread == 0 demo tsukau
{
	compactdistarrthread_arg_t *targ = (compactdistarrthread_arg_t *)arg;
	int njob = targ->njob;
	int para = targ->para;
	int im = targ->im;
	int jm = targ->jm;
//	int thread_no = targ->thread_no;
	int *nlen = targ->nlen;
	int **pointt = targ->pointt;
	int *table1 = targ->table1;
	int *table2 = targ->table2;
	int *tselfscore = targ->tselfscore;
 	int *joblist = targ->joblist;
 	int *posshared = targ->posshared;
	double *mindist = targ->mindist;
	int *nearest = targ->nearest;
//	double **partmtx = targ->partmtx;
	double *newarr = targ->newarr;
	int i, posinjoblist, n;

	double tmpdist1;
	double tmpdist2;
	double tmpdouble;

//			for( acpti=ac; acpti!=NULL; acpti=acpti->next )

	while( 1 )
	{
#ifdef enablemultithread
		if( para ) pthread_mutex_lock( targ->mutex );
#endif 
		if( *posshared >= njob ) // block no toki >=
		{
#ifdef enablemultithread
			if( para ) pthread_mutex_unlock( targ->mutex );
#endif
			commonsextet_p( NULL, NULL );
			return( NULL );
		}
		posinjoblist = *posshared;
		*posshared += BLOCKSIZE;
#ifdef enablemultithread
		if( para ) pthread_mutex_unlock( targ->mutex );
#endif

		for( n=0; n<BLOCKSIZE&&posinjoblist<njob; n++ )
		{
			i = joblist[posinjoblist++];

			if( i == im ) continue;
			if( i == jm ) continue;
	
//			if( partmtx[im] )
//				tmpdist1 = partmtx[im][i];
//			else if( partmtx[i] )
//				tmpdist1 = partmtx[i][im];
//			else
				tmpdist1 = distcompact( nlen[im], nlen[i], table1, pointt[i], tselfscore[im], tselfscore[i] );
					
//			if( partmtx[jm] )
//				tmpdist2 = partmtx[jm][i];
//			else if( partmtx[i] )
//				tmpdist2 = partmtx[i][jm];
//			else
				tmpdist2 = distcompact( nlen[jm], nlen[i], table2, pointt[i], tselfscore[jm], tselfscore[i] );
	
//			if( seq )
//			{
//				tmpdist1 = distcompact_msa( seq[im], seq[i], skiptable[im], skiptable[i], tselfscore[im], tselfscore[i] );
//				tmpdist2 = distcompact_msa( seq[jm], seq[i], skiptable[jm], skiptable[i], tselfscore[jm], tselfscore[i] );
//			}
//			else
//			{
//				tmpdist1 = distcompact( nlen[im], nlen[i], table1, pointt[i], tselfscore[im], tselfscore[i] );
//				tmpdist2 = distcompact( nlen[jm], nlen[i], table2, pointt[i], tselfscore[jm], tselfscore[i] );
//			}
			tmpdouble = cluster_mix_double( tmpdist1, tmpdist2 );
			newarr[i] = tmpdouble;
	
//			if( partmtx[i] ) partmtx[i][im] = partmtx[i][jm] = newarr[i];
	
			if( tmpdouble < mindist[i]  )
			{
				mindist[i] = tmpdouble;
				nearest[i] = im;
			}
	
//			if( tmpdouble < mindist[im]  ) // koko deha muri
//			{
//				mindist[im] = tmpdouble;
//				nearest[im] = i;
//			}
	
			if( nearest[i] == jm )
			{
				nearest[i] = im;
			}
		}
	}
}

static void *kmerdistarrthreadjoblist( void *arg ) // enablemultithread == 0 demo tsukau
{
	compactdistarrthread_arg_t *targ = (compactdistarrthread_arg_t *)arg;
	int njob = targ->njob;
	int para = targ->para;
	int im = targ->im;
	int jm = targ->jm;
//	int thread_no = targ->thread_no;
	int *nlen = targ->nlen;
	int **pointt = targ->pointt;
	int *table1 = targ->table1;
	int *table2 = targ->table2;
	int *tselfscore = targ->tselfscore;
 	int *joblist = targ->joblist;
 	int *posshared = targ->posshared;
	double *mindist = targ->mindist;
	int *nearest = targ->nearest;
	double **partmtx = targ->partmtx;
	double *newarr = targ->newarr;
	int i, posinjoblist, n;

	double tmpdist1;
	double tmpdist2;
	double tmpdouble;

//			for( acpti=ac; acpti!=NULL; acpti=acpti->next )

	while( 1 )
	{
#ifdef enablemultithread
		if( para ) pthread_mutex_lock( targ->mutex );
#endif 
		if( *posshared >= njob ) // block no toki >=
		{
#ifdef enablemultithread
			if( para ) pthread_mutex_unlock( targ->mutex );
#endif
			commonsextet_p( NULL, NULL );
			return( NULL );
		}
		posinjoblist = *posshared;
		*posshared += BLOCKSIZE;
#ifdef enablemultithread
		if( para ) pthread_mutex_unlock( targ->mutex );
#endif

		for( n=0; n<BLOCKSIZE&&posinjoblist<njob; n++ )
		{
			i = joblist[posinjoblist++];

			if( i == im ) continue;
			if( i == jm ) continue;
	
			if( partmtx[im] )
				tmpdist1 = partmtx[im][i];
			else if( partmtx[i] )
				tmpdist1 = partmtx[i][im];
			else
				tmpdist1 = distcompact( nlen[im], nlen[i], table1, pointt[i], tselfscore[im], tselfscore[i] );
					
			if( partmtx[jm] )
				tmpdist2 = partmtx[jm][i];
			else if( partmtx[i] )
				tmpdist2 = partmtx[i][jm];
			else
				tmpdist2 = distcompact( nlen[jm], nlen[i], table2, pointt[i], tselfscore[jm], tselfscore[i] );
	
//			if( seq )
//			{
//				tmpdist1 = distcompact_msa( seq[im], seq[i], skiptable[im], skiptable[i], tselfscore[im], tselfscore[i] );
//				tmpdist2 = distcompact_msa( seq[jm], seq[i], skiptable[jm], skiptable[i], tselfscore[jm], tselfscore[i] );
//			}
//			else
//			{
//				tmpdist1 = distcompact( nlen[im], nlen[i], table1, pointt[i], tselfscore[im], tselfscore[i] );
//				tmpdist2 = distcompact( nlen[jm], nlen[i], table2, pointt[i], tselfscore[jm], tselfscore[i] );
//			}
			tmpdouble = cluster_mix_double( tmpdist1, tmpdist2 );
			newarr[i] = tmpdouble;
	
			if( partmtx[i] ) partmtx[i][im] = partmtx[i][jm] = newarr[i];
	
			if( tmpdouble < mindist[i]  )
			{
				mindist[i] = tmpdouble;
				nearest[i] = im;
			}
	
//			if( tmpdouble < mindist[im]  ) // koko deha muri
//			{
//				mindist[im] = tmpdouble;
//				nearest[im] = i;
//			}
	
			if( nearest[i] == jm )
			{
				nearest[i] = im;
			}
		}
	}
}

static void *verycompactmsadistarrthreadjoblist( void *arg ) // enablemultithread == 0 demo tsukau
{
	compactdistarrthread_arg_t *targ = (compactdistarrthread_arg_t *)arg;
	int njob = targ->njob;
	int para = targ->para;
	int im = targ->im;
	int jm = targ->jm;
//	int thread_no = targ->thread_no;
	int *tselfscore = targ->tselfscore;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
 	int *joblist = targ->joblist;
 	int *posshared = targ->posshared;
	double *mindist = targ->mindist;
	int *nearest = targ->nearest;
//	double **partmtx = targ->partmtx;
	double *newarr = targ->newarr;
	int i, posinjoblist, n;

	double tmpdist1;
	double tmpdist2;
	double tmpdouble;

//			for( acpti=ac; acpti!=NULL; acpti=acpti->next )


	while( 1 )
	{
#ifdef enablemultithread
		if( para ) pthread_mutex_lock( targ->mutex );
#endif 
		if( *posshared >= njob ) // block no toki >=
		{
#ifdef enablemultithread
			if( para ) pthread_mutex_unlock( targ->mutex );
#endif
			commonsextet_p( NULL, NULL );
			return( NULL );
		}
		posinjoblist = *posshared;
		*posshared += BLOCKSIZE;
#ifdef enablemultithread
		if( para ) pthread_mutex_unlock( targ->mutex );
#endif

		for( n=0; n<BLOCKSIZE&&posinjoblist<njob; n++ )
		{
			i = joblist[posinjoblist++];

			if( i == im ) continue;
			if( i == jm ) continue;

//			if( partmtx[im] )
//				tmpdist1 = partmtx[im][i];
//			else if( partmtx[i] )
//				tmpdist1 = partmtx[i][im];
//			else
				tmpdist1 = distcompact_msa( seq[im], seq[i], skiptable[im], skiptable[i], tselfscore[im], tselfscore[i] );
					
//			if( partmtx[jm] )
//				tmpdist2 = partmtx[jm][i];
//			else if( partmtx[i] )
//				tmpdist2 = partmtx[i][jm];
//			else
				tmpdist2 = distcompact_msa( seq[jm], seq[i], skiptable[jm], skiptable[i], tselfscore[jm], tselfscore[i] );
	
			tmpdouble = cluster_mix_double( tmpdist1, tmpdist2 );
			newarr[i] = tmpdouble;
	
//			if( partmtx[i] ) partmtx[i][im] = partmtx[i][jm] = newarr[i];
	
			if( tmpdouble < mindist[i]  )
			{
				mindist[i] = tmpdouble;
				nearest[i] = im;
			}
	
//			if( tmpdouble < mindist[im]  ) // koko deha muri
//			{
//				mindist[im] = tmpdouble;
//				nearest[im] = i;
//			}
	
			if( nearest[i] == jm )
			{
				nearest[i] = im;
			}
		}
	}
}

static void *msadistarrthreadjoblist( void *arg ) // enablemultithread == 0 demo tsukau
{
	compactdistarrthread_arg_t *targ = (compactdistarrthread_arg_t *)arg;
	int njob = targ->njob;
	int para = targ->para;
	int im = targ->im;
	int jm = targ->jm;
//	int thread_no = targ->thread_no;
	int *tselfscore = targ->tselfscore;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
 	int *joblist = targ->joblist;
 	int *posshared = targ->posshared;
	double *mindist = targ->mindist;
	int *nearest = targ->nearest;
	double **partmtx = targ->partmtx;
	double *newarr = targ->newarr;
	int i, posinjoblist, n;

	double tmpdist1;
	double tmpdist2;
	double tmpdouble;

//			for( acpti=ac; acpti!=NULL; acpti=acpti->next )


	while( 1 )
	{
#ifdef enablemultithread
		if( para ) pthread_mutex_lock( targ->mutex );
#endif 
		if( *posshared >= njob ) // block no toki >=
		{
#ifdef enablemultithread
			if( para ) pthread_mutex_unlock( targ->mutex );
#endif
			commonsextet_p( NULL, NULL );
			return( NULL );
		}
		posinjoblist = *posshared;
		*posshared += BLOCKSIZE;
#ifdef enablemultithread
		if( para ) pthread_mutex_unlock( targ->mutex );
#endif

		for( n=0; n<BLOCKSIZE&&posinjoblist<njob; n++ )
		{
			i = joblist[posinjoblist++];

			if( i == im ) continue;
			if( i == jm ) continue;

			if( partmtx[im] )
				tmpdist1 = partmtx[im][i];
			else if( partmtx[i] )
				tmpdist1 = partmtx[i][im];
			else
				tmpdist1 = distcompact_msa( seq[im], seq[i], skiptable[im], skiptable[i], tselfscore[im], tselfscore[i] );
					
			if( partmtx[jm] )
				tmpdist2 = partmtx[jm][i];
			else if( partmtx[i] )
				tmpdist2 = partmtx[i][jm];
			else
				tmpdist2 = distcompact_msa( seq[jm], seq[i], skiptable[jm], skiptable[i], tselfscore[jm], tselfscore[i] );
	
			tmpdouble = cluster_mix_double( tmpdist1, tmpdist2 );
			newarr[i] = tmpdouble;
	
			if( partmtx[i] ) partmtx[i][im] = partmtx[i][jm] = newarr[i];
	
			if( tmpdouble < mindist[i]  )
			{
				mindist[i] = tmpdouble;
				nearest[i] = im;
			}
	
//			if( tmpdouble < mindist[im]  ) // koko deha muri
//			{
//				mindist[im] = tmpdouble;
//				nearest[im] = i;
//			}
	
			if( nearest[i] == jm )
			{
				nearest[i] = im;
			}
		}
	}
}

void compacttree_memsaveselectable( int nseq, double **partmtx, int *nearest, double *mindist, int **pointt, int *tselfscore, char **seq, int **skiptable, int ***topol, double **len, char **name, int *nlen, Treedep *dep, int treeout, int howcompact, int memsave )
{
	int i, j, k;
//	int miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
//	double tmpdouble;
//	double eff1, eff0;
	double *tmptmplen = NULL; //static?
	int *hist = NULL; //static?
	Bchain *ac = NULL; //static?
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti;
	int *pt1, *pt2, *pt11, *pt22;
	int *nmemar; //static?
	int nmemim, nmemjm;
	double minscore;
	char **tree; //static?
	char *treetmp; //static?
	char *nametmp, *nameptr, *tmpptr; //static?
	FILE *fp;
	double (*clusterfuncpt[1])(double,double);
	char namec;
	int *singlettable1 = NULL;
	int *singlettable2 = NULL;
	double *newarr;
	void *(*distarrfunc)( void * );
	void *(*resetnearestfunc)( void * );
	int numfilled;
	int nthreadtree;
	compactdistarrthread_arg_t *distarrarg;
	resetnearestthread_arg_t *resetarg;
	int *joblist, nactive, posshared;
	double *result;


	sueff1 = 1 - (double)sueff_global;
	sueff05 = (double)sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_double;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}

	if( howcompact == 2 )
	{
		if( seq )
		{
//			distarrfunc = verycompactmsadistarrthread;
			distarrfunc = verycompactmsadistarrthreadjoblist;
			resetnearestfunc = NULL;
		}
		else
		{
//			distarrfunc = verycompactkmerdistarrthread;
			distarrfunc = verycompactkmerdistarrthreadjoblist;
			resetnearestfunc = NULL;
		}
	}
	else
	{
		if( seq )
		{
			distarrfunc = msadistarrthreadjoblist;
			resetnearestfunc = msaresetnearestthread;
		}	
		else
		{
			distarrfunc = kmerdistarrthreadjoblist;
			resetnearestfunc = kmerresetnearestthread;
		}
	}
	distarrarg = calloc( MAX( nthread, 1 ), sizeof( compactdistarrthread_arg_t ) );
	resetarg = calloc( MAX( nthread, 1 ), sizeof( resetnearestthread_arg_t ) );
	joblist = calloc( njob, sizeof( int ) );
	if( howcompact != 2 ) result = calloc( njob, sizeof( double ) );
	else result = NULL;

	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateFloatVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		if( treeout )
		{
			treetmp = NULL; // kentou 2013/06/12
			nametmp = AllocateCharVec( 1000 ); // nagasugi
			tree = AllocateCharMtx( njob, 0 );
		}
	}

	
	if( treeout )
	{
	    for( i=0; i<nseq; i++ )
		{
			for( j=0; j<999; j++ ) nametmp[j] = 0;
			for( j=0; j<999; j++ ) 
			{
				namec = name[i][j];
				if( namec == 0 )
					break;
				else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
					nametmp[j] = namec;
				else
					nametmp[j] = '_';
			}
			nametmp[j] = 0;
//			sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
			if( outnumber )
				nameptr = strstr( nametmp, "_numo_e" ) + 8;
			else
				nameptr = nametmp + 1;
	
			if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame
	
			tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
			if( tree[i] == NULL )
			{
				reporterr(       "Cannot allocate tree!\n" );
				exit( 1 );
			}
			sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
		}
	}

	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

//	for( i=0; i<nseq; i++ ) setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	for( i=0,numfilled=0; i<nseq; i++ ) if( partmtx[i] ) numfilled++;
	reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 100 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );

//		for( i=0,j=0; i<nseq; i++ ) if( partmtx[i] ) j++;
//		if( k% 100 == 0 ) reporterr( "numfilled=%d, filledinpartmtx=%d, numempty=%d\n", numfilled, j, nseq-k-numfilled );

		minscore = 999.9;
		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
		{
			i = acpti->pos;
//			printf(       "k=%d i=%d, mindist[i]=%f\n", k, i, mindist[i] );
			if( mindist[i] < minscore ) // muscle
			{
				im = i;
				minscore = mindist[i];
			}
		}
//		printf(       "minscore=%f\n", minscore );
		jm = nearest[im];
//		printf(       "im=%d\n", im );
//		printf(       "jm=%d\n", jm );


		if( jm < im ) 
		{
			j=jm; jm=im; im=j;
		}

		if( partmtx[im] == NULL && howcompact != 2 ) numfilled++;
		if( partmtx[jm] != NULL ) numfilled--;

		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];
		if( memsave )
			intpt = topol[k][0] = (int *)realloc( topol[k][0], ( 2 ) * sizeof( int ) ); // memsave
		else
			intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) ); // memsave
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
//				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
//				pt22 = pt2;
			}
			if( memsave )
			{
				*intpt++ = *pt11;
				*intpt = -1;
			}
			else
			{
				for( intpt2=pt11; *intpt2!=-1; )
					*intpt++ = *intpt2++;
				for( intpt2=pt22; *intpt2!=-1; )
					*intpt++ = *intpt2++;
				*intpt = -1;
			}
		}

		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;
		nmemjm = nmemar[jm];
		if( memsave )
			intpt = topol[k][1] = (int *)realloc( topol[k][1], ( 2 ) * sizeof( int ) ); // memsave
		else
			intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) ); // memsave
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
//				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
//				pt22 = pt2;
			}
			if( memsave )
			{
				*intpt++ = *pt11;
				*intpt = -1;
			}
			else
			{
				for( intpt2=pt11; *intpt2!=-1; )
					*intpt++ = *intpt2++;
				for( intpt2=pt22; *intpt2!=-1; )
					*intpt++ = *intpt2++;
				*intpt = -1;
			}
		}

		minscore *= 0.5;

//		printf( "minscore = %f, tmptmplen[im] = %f, tmptmplen[jm] = %f\n", minscore, tmptmplen[im], tmptmplen[jm] );

		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );

		if( dep ) dep[k].distfromtip = minscore;
//		reporterr(       "\n##### dep[%d].distfromtip = %f\n", k, minscore );

		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;
		mindist[im] = 999.9;

		if( pointt ) // kmer
		{
			singlettable1 = (int *)calloc( tsize, sizeof( int ) );
			singlettable2 = (int *)calloc( tsize, sizeof( int ) );
			makecompositiontable_global( singlettable1, pointt[im] );
			makecompositiontable_global( singlettable2, pointt[jm] );
		}

		newarr = calloc( nseq, sizeof( double ) );

//		nthreadtree = MAX( 1, nthread );
		nthreadtree = nthread;


		for( acpti=ac,nactive=0; acpti!=NULL; acpti=acpti->next ) joblist[nactive++] = acpti->pos; // sukoshi muda...


#ifdef enablemultithread
		if( nthreadtree > 0 )
		{
			compactdistarrthread_arg_t *targ;
			pthread_t *handle;
			pthread_mutex_t mutex;
		
			posshared = 0;
//			targ = calloc( nthreadtree, sizeof( compactdistarrthread_arg_t ) );
			targ = distarrarg;
			handle = calloc( nthreadtree, sizeof( pthread_t ) );
			pthread_mutex_init( &mutex, NULL );
		
			if( k % 100 == 0 ) reporterr( " (%d threads, nactive=%d, nfilled=%d)     \r", nthreadtree, nactive, numfilled );
			for( i=0; i<nthreadtree; i++ )
			{
				targ[i].para = 1;
				targ[i].njob = nactive;
//				targ[i].thread_no = i;
				targ[i].im = im;
				targ[i].jm = jm;
				targ[i].tselfscore = tselfscore;
				targ[i].nlen = nlen;
				targ[i].seq = seq;
				targ[i].skiptable = skiptable;
				targ[i].pointt = pointt;
				targ[i].table1 = singlettable1;
				targ[i].table2 = singlettable2;
				targ[i].joblist = joblist;
				targ[i].posshared = &posshared;
				targ[i].mindist = mindist;
				targ[i].nearest = nearest;
				targ[i].newarr = newarr;
				targ[i].partmtx = partmtx;
				targ[i].mutex = &mutex;

				pthread_create( handle+i, NULL, distarrfunc, (void *)(targ+i) );
			}
		
			for( j=0; j<nthreadtree; j++ ) pthread_join( handle[j], NULL );
			pthread_mutex_destroy( &mutex );
			free( handle );
//			free( targ );
	
#if 0
			for( acpti=ac; acpti!=NULL; acpti=acpti->next ) // antei sei no tame
			{
				i = acpti->pos;
				if( i != im && i != jm )
				{
//					if( partmtx[i] ) partmtx[i][im] = partmtx[i][jm] = newarr[i]; // heiretsu demo ii.
//					if( newarr[i] < mindist[i]  )
//					{
//						mindist[i] = newarr[i];
//						nearest[i] = im;
//					}
					if( newarr[i] < mindist[im]  )
					{
						mindist[im] = newarr[i];
						nearest[im] = i;
					}
//					if( nearest[i] == jm )
//					{
//						nearest[i] = im;
//					}
				}
			}
#endif
		}
		else
#endif
		{
			if( k % 100 == 0 ) reporterr( " (serial, nactive=%d, nfilled=%d)             \r", nactive, numfilled );
			compactdistarrthread_arg_t *targ;
		
			posshared = 0;
//			targ = calloc( 1, sizeof( compactdistarrthread_arg_t ) );
			targ = distarrarg;
		
			for( i=0; i<1; i++ )
			{
				targ[i].para = 0;
				targ[i].njob = nactive;
//				targ[i].thread_no = i;
				targ[i].im = im;
				targ[i].jm = jm;
				targ[i].tselfscore = tselfscore;
				targ[i].nlen = nlen;
				targ[i].seq = seq;
				targ[i].skiptable = skiptable;
				targ[i].pointt = pointt;
				targ[i].table1 = singlettable1;
				targ[i].table2 = singlettable2;
				targ[i].joblist = joblist;
				targ[i].posshared = &posshared;
				targ[i].mindist = mindist;
				targ[i].nearest = nearest;
				targ[i].newarr = newarr;
				targ[i].partmtx = partmtx;

				distarrfunc( targ+i );
//				pthread_create( handle, NULL, distarrfunc, (void *)(targ) );
			}

//			free( targ );
	
		}

		for( acpti=ac; acpti!=NULL; acpti=acpti->next ) // antei sei no tame
		{
			i = acpti->pos;
			if( i != im && i != jm )
			{
//				if( partmtx[i] ) partmtx[i][im] = partmtx[i][jm] = newarr[i]; // heiretsu demo ii.
//				if( newarr[i] < mindist[i]  )
//				{
//					mindist[i] = newarr[i];
//					nearest[i] = im;
//				}
				if( newarr[i] < mindist[im]  )
				{
					mindist[im] = newarr[i];
					nearest[im] = i;
				}
//				if( nearest[i] == jm )
//				{
//					nearest[i] = im;
//				}
			}
		}

//		printf( "im=%d, jm=%d\n", im, jm );
#if 0
		printf( "matrix = \n" );
		for( i=0; i<njob; i++ )
		{
			if( partmtx[i] ) for( j=0; j<njob; j++ ) printf( "%f ", partmtx[i][j] );
			else printf( "nai" );
			printf( "\n" );
			
		}
#endif
//		if( k%500 == 0 )
//		{
//			reporterr( "at step %d,", k );
//			use_getrusage();
//		}

		if( partmtx[im] ) free( partmtx[im] ); partmtx[im] = NULL;
		if( partmtx[jm] ) free( partmtx[jm] ); partmtx[jm] = NULL;
		if( howcompact == 2 )
		{
			free( newarr );
			newarr = NULL;
		}
		else
		{
			partmtx[im] = newarr;
		}


		if( pointt )
		{
			free( singlettable1 );
			free( singlettable2 );
			singlettable1 = NULL;
			singlettable2 = NULL;
		}
	
		if( treeout )
		{
			treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
			if( !treetmp )
			{
				reporterr(       "Cannot allocate treetmp\n" );
				exit( 1 );
			}
			sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
			free( tree[im] );
			free( tree[jm] );
			tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
			tree[jm] = NULL;
			if( tree[im] == NULL )
			{
				reporterr(       "Cannot reallocate tree!\n" );
				exit( 1 );
			}
			strcpy( tree[im], treetmp );
		}

		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;

#if 0 // muscle seems to miss this.
//		int nwork = 0;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
//			printf( "reset nearest? i=%d, k=%d, nearest[i]=%d, im=%d, mindist=%f\n", i, k, nearest[i], im, mindist[i] );
			if( nearest[i] == im ) 
			{
//				printf( "reset nearest, i=%d, k=%d\n", i, k );
				if( partmtx[im][i] > mindist[i] )
				{
//					nwork++;
//					printf( "go\n" );
					if( pointt ) // kmer
					{
						singlettable1 = (int *)calloc( tsize, sizeof( int ) );
						makecompositiontable_global( singlettable1, pointt[i] );
					}
					resetnearest( nseq, ac, partmtx, mindist+i, nearest+i, i, seq, skiptable, tselfscore, pointt, nlen, singlettable1 );
					if( pointt ) free( singlettable1 ); singlettable1 = NULL;// kmer
					if( pointt ) commonsextet_p( NULL, NULL );
				}
			}
		}
//		reporterr( "nwork = %d\n", nwork );
#else

		if( howcompact == 2 ) continue;

#if 0
		if( 0 && nthreadtree > 0  )
		{
			resetnearestthread_arg_t *targ;
			pthread_t *handle;
			pthread_mutex_t mutex;
			Bchain *acshared;
		
			acshared = ac;
//			targ = calloc( nthreadtree, sizeof( resetnearestthread_arg_t ) );
			targ = resetarg;
			handle = calloc( nthreadtree, sizeof( pthread_t ) );
			pthread_mutex_init( &mutex, NULL );
		
			for( i=0; i<nthreadtree; i++ )
			{
				targ[i].para = 1;
				targ[i].nseq = nseq;
				targ[i].im = im;
				targ[i].partmtx = partmtx;
				targ[i].mindist = mindist;
				targ[i].nearest = nearest;
				targ[i].seq = seq;
				targ[i].skiptable = skiptable;
				targ[i].tselfscore = tselfscore;
				targ[i].pointt = pointt;
				targ[i].nlen = nlen;
				targ[i].acpt = &acshared;
				targ[i].ac = ac;
				targ[i].mutex = &mutex;

				pthread_create( handle+i, NULL, resetnearestfunc, (void *)(targ+i) );
			}
		
			for( j=0; j<nthreadtree; j++ ) pthread_join( handle[j], NULL );
			pthread_mutex_destroy( &mutex );
			free( handle );
//			free( targ );
		}
		else
#endif
		{
			Bchain *acshared;
			acshared = ac;
			resetnearestthread_arg_t *targ;
//			targ = calloc( 1, sizeof( resetnearestthread_arg_t ) );
			targ = resetarg;
			{
				targ[0].para = 0;
				targ[0].nseq = nseq;
				targ[0].im = im;
				targ[0].partmtx = partmtx;
				targ[0].mindist = mindist;
				targ[0].nearest = nearest;
				targ[0].seq = seq;
				targ[0].skiptable = skiptable;
				targ[0].tselfscore = tselfscore;
				targ[0].pointt = pointt;
				targ[0].nlen = nlen;
				targ[0].result = result;
				targ[0].joblist = joblist;
				targ[0].acpt = &acshared;
				targ[0].ac = ac;

				resetnearestfunc( targ );
			}
//			free( targ );
		}
#endif


#if 0
        printf(       "\nooSTEP-%03d:\n", k+1 );
		printf(       "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) printf(       " %03d", topol[k][0][i]+1 );
        printf(       "\n" );
		printf(       "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) printf(       " %03d", topol[k][1][i]+1 );
        printf(       "\n" );
#endif
    }
	if( treeout )
	{
		fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
		fclose( fp );
	}

	for( im=0; im<nseq; im++ ) // im wo ugokasu hituyouha nai.
	{
		if( partmtx[im] ) free( partmtx[im] ); partmtx[im] = NULL;
	}
//	if( partmtx ) free( partmtx ); partmtx = NULL; // oya ga free
	if( treeout )
	{
		free( tree[0] );
		free( tree );
		free( treetmp );
		free( nametmp );
	}
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	if( singlettable1 ) free( singlettable1 );
	if( singlettable2 ) free( singlettable2 );
	free( distarrarg );
	free( resetarg );
	free( joblist );
	if( result ) free( result );
}

void fixed_musclesupg_double_realloc_nobk_halfmtx_treeout_memsave( int nseq, double **eff, int ***topol, double **len, char **name, int *nlen, Treedep *dep, int efffree )
{

	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt;
	double tmpdouble;
	double eff1, eff0;
	double *tmptmplen = NULL; //static?
	int *hist = NULL; //static?
	Bchain *ac = NULL; //static?
	int im = 1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti;
	int *pt1, *pt2, *pt11;
	int *nmemar; //static?
	int nmemim, nmemjm;
	double minscore;
	int *nearest = NULL; // by D.Mathog, a guess
	double *mindisfrom = NULL; // by D.Mathog, a guess
	char **tree; //static?
	char *treetmp; //static?
	char *nametmp, *nameptr, *tmpptr; //static?
	FILE *fp;
	double (*clusterfuncpt[1])(double,double);
	char namec;


	sueff1 = 1 - (double)sueff_global;
	sueff05 = (double)sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_double;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_double;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_double;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}

	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateFloatVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateFloatVec( njob );
		nearest = AllocateIntVec( njob );
//		treetmp = AllocateCharVec( njob * ( B + 100 ) ); // nagasugi?
		treetmp = NULL; // kentou 2013/06/12
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( njob, njob*600 );
		tree = AllocateCharMtx( njob, 0 );
	}

	
    for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

	for( i=0; i<nseq; i++ ) setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle


	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );

		minscore = 999.9;
		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
		{
			i = acpti->pos;
//			printf(       "k=%d i=%d, mindist[i]=%f\n", k, i, mindisfrom[i] );
			if( mindisfrom[i] < minscore ) // muscle
			{
				im = i;
				minscore = mindisfrom[i];
			}
		}

//		printf(       "minscore=%f\n", minscore );
		jm = nearest[im];
//		printf(       "im=%d\n", im );
//		printf(       "jm=%d\n", jm );
		if( jm < im ) 
		{
			j=jm; jm=im; im=j;
		}


		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];
		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( 2 ) * sizeof( int ) ); // memsave
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
//				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
//				pt22 = pt2;
			}
#if 1 // memsave
			*intpt++ = *pt11;
			*intpt = -1;
#else
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
#endif
		}

		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;
		nmemjm = nmemar[jm];
		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( 2 ) * sizeof( int ) ); // memsave
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
//				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
//				pt22 = pt2;
			}
#if 1 // memsave
			*intpt++ = *pt11;
			*intpt = -1;
#else
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
#endif
		}

		minscore *= 0.5;

//		printf( "minscore = %f, tmptmplen[im] = %f, tmptmplen[jm] = %f\n", minscore, tmptmplen[im], tmptmplen[jm] );
		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );

		if( dep ) dep[k].distfromtip = minscore;
//		reporterr(       "\n##### dep[%d].distfromtip = %f\n", k, minscore );

		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

		mindisfrom[im] = 999.9;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim-miniim];
				eff1 = eff[minijm][maxijm-minijm];
#if 0
                		tmpdouble = eff[miniim][maxiim-miniim] =
				MIN( eff0, eff1 ) * sueff1 + ( eff0 + eff1 ) * sueff05; 
#else
                tmpdouble = eff[miniim][maxiim-miniim] =
				(clusterfuncpt[0])( eff0, eff1 );
//				printf( "tmpdouble=%f, eff0=%f, eff1=%f\n", tmpdouble, eff0, eff1 );
#endif
				if( tmpdouble < mindisfrom[i]  )
				{
					mindisfrom[i] = tmpdouble;
					nearest[i] = im;
				}
				if( tmpdouble < mindisfrom[im]  )
				{
					mindisfrom[im] = tmpdouble;
					nearest[im] = i;
				}
				if( nearest[i] == jm )
				{
					nearest[i] = im;
				}
            }
        }
//		printf( "im=%d, jm=%d\n", im, jm );
#if 0
		printf( "matrix = \n" );
		for( i=0; i<njob; i++ )
		{
			for( j=0; j<njob; j++ ) 
			{
				if( i>j )
				{
					minijm=j;
					maxijm=i;
				}
				else
				{
					minijm=i;
					maxijm=j;
				}
				printf( "%f ", eff[minijm][maxijm-minijm] );
			}
			printf( "\n" );
		}
#endif

		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );

		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
		if( efffree )
		{
			free( (void *)eff[jm] ); eff[jm] = NULL; // Ato de fukkatsu
		}

#if 1 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
//			printf( "reset nearest? i=%d, k=%d, nearest[i]=%d, im=%d, mindist=%f\n", i, k, nearest[i], im, mindisfrom[i] );
			if( nearest[i] == im ) 
			{
//				printf( "reset nearest, i=%d, k=%d\n", i, k );
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
				}
				if( eff[miniim][maxiim-miniim] > mindisfrom[i] )
				{
//					printf( "go\n" );
					setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i );
				}
			}
		}
#else
		reporterr( "CHUUI!\n" );
#endif


#if 0
        printf(       "\nooSTEP-%03d:\n", k+1 );
		printf(       "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) printf(       " %03d", topol[k][0][i]+1 );
        printf(       "\n" );
		printf(       "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) printf(       " %03d", topol[k][1][i]+1 );
        printf(       "\n" );
#endif
    }
	fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
	fclose( fp );

	free( tree[0] );
	free( tree );
	free( treetmp );
	free( nametmp );
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	free( mindisfrom );
	free( nearest );
}

void fixed_musclesupg_double_realloc_nobk_halfmtx_treeout( int nseq, double **eff, int ***topol, double **len, char **name, int *nlen, Treedep *dep, int efffree )
{
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	double tmpdouble;
	double eff1, eff0;
	double *tmptmplen = NULL; //static?
	int *hist = NULL; //static?
	Bchain *ac = NULL; //static?
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti;
	int *pt1, *pt2, *pt11, *pt22;
	int *nmemar; //static?
	int nmemim, nmemjm;
	double minscore;
	int *nearest = NULL; // by D.Mathog, a guess
	double *mindisfrom = NULL; // by D.Mathog, a guess
	char **tree; //static?
	char *treetmp; //static?
	char *nametmp, *nameptr, *tmpptr; //static?
	FILE *fp;
	double (*clusterfuncpt[1])(double,double);
	char namec;


	sueff1 = 1 - (double)sueff_global;
	sueff05 = (double)sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_double;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_double;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_double;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}

	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateFloatVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateFloatVec( njob );
		nearest = AllocateIntVec( njob );
//		treetmp = AllocateCharVec( njob * ( B + 100 ) ); // nagasugi?
		treetmp = NULL; // kentou 2013/06/12
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( njob, njob*600 );
		tree = AllocateCharMtx( njob, 0 );
	}

	
    for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

	for( i=0; i<nseq; i++ ) setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );

		minscore = 999.9;
		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
		{
			i = acpti->pos;
//			reporterr(       "k=%d i=%d\n", k, i );
			if( mindisfrom[i] < minscore ) // muscle
			{
				im = i;
				minscore = mindisfrom[i];
			}
		}
		jm = nearest[im];
		if( jm < im ) 
		{
			j=jm; jm=im; im=j;
		}


		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];
		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;
		nmemjm = nmemar[jm];
		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );

		if( dep ) dep[k].distfromtip = minscore;
//		reporterr(       "\n##### dep[%d].distfromtip = %f\n", k, minscore );

		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

		mindisfrom[im] = 999.9;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim-miniim];
				eff1 = eff[minijm][maxijm-minijm];
#if 0
                		tmpdouble = eff[miniim][maxiim-miniim] =
				MIN( eff0, eff1 ) * sueff1 + ( eff0 + eff1 ) * sueff05; 
#else
                tmpdouble = eff[miniim][maxiim-miniim] =
				(clusterfuncpt[0])( eff0, eff1 );


#endif
				if( tmpdouble < mindisfrom[i]  )
				{
					mindisfrom[i] = tmpdouble;
					nearest[i] = im;
				}
				if( tmpdouble < mindisfrom[im]  )
				{
					mindisfrom[im] = tmpdouble;
					nearest[im] = i;
				}
				if( nearest[i] == jm )
				{
					nearest[i] = im;
				}
            }
        }

		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );

		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
		if( efffree )
		{
			free( (void *)eff[jm] ); eff[jm] = NULL;
		}

#if 1 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
				}
				if( eff[miniim][maxiim-miniim] > mindisfrom[i] )
					setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#else
		reporterr( "chuui!\n" );
#endif


#if 0
        printf(       "\nooSTEP-%03d:\n", k+1 );
		printf(       "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) printf(       " %03d", topol[k][0][i]+1 );
        printf(       "\n" );
		printf(       "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) printf(       " %03d", topol[k][1][i]+1 );
        printf(       "\n" );
#endif
    }
	fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
	fclose( fp );

	free( tree[0] );
	free( tree );
	free( treetmp );
	free( nametmp );
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	free( mindisfrom );
	free( nearest );
}

void fixed_musclesupg_double_treeout( int nseq, double **eff, int ***topol, double **len, char **name )
{
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	double tmpdouble;
	double eff1, eff0;
	static double *tmptmplen = NULL;
	static int *hist = NULL;
	static Bchain *ac = NULL;
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti;
	int *pt1, *pt2, *pt11, *pt22;
	static int *nmemar;
	int nmemim, nmemjm;
	double minscore;
	int *nearest = NULL; // by D.Mathog, a guess
	double *mindisfrom = NULL; // by D.Mathog, a guess
	static char **tree;
	static char *treetmp;
	static char *nametmp, *nameptr, *tmpptr;
	FILE *fp;
	double (*clusterfuncpt[1])(double,double);
	char namec;


	sueff1 = 1.0 - sueff_global;
	sueff05 = sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_double;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_double;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_double;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}





#if 0
	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateDoubleVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateDoubleVec( njob );
		nearest = AllocateIntVec( njob );
		treetmp = AllocateCharVec( njob*150 );
		nametmp = AllocateCharVec( 91 );
		tree = AllocateCharMtx( njob, njob*150 );
	}
    for( i=0; i<nseq; i++ )
	{
		for( j=0; j<90; j++ ) nametmp[j] = 0;
		for( j=0; j<90; j++ ) 
		{
			if( name[i][j] == 0 )
				break;
			else if( isalnum( name[i][j] ) )
				nametmp[j] = name[i][j];
			else
				nametmp[j] = '_';
		}
		nametmp[90] = 0;
//		sprintf( tree[i], "%d_%.60s", i+1, nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		sprintf( tree[i], "\n%d_%.60s\n", i+1, nameptr );
	}

#else

	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateDoubleVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateDoubleVec( njob );
		nearest = AllocateIntVec( njob );
//		treetmp = AllocateCharVec( njob * ( B + 100 ) ); // nagasugi?
		treetmp = NULL; // kentou 2013/06/12
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( njob, njob*600 );
		tree = AllocateCharMtx( njob, 0 );
	}

	
    for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}

#endif








	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

	for( i=0; i<nseq; i++ ) setnearest_double_fullmtx( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );

		minscore = 999.9;
		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
		{
			i = acpti->pos;
//			reporterr(       "k=%d i=%d\n", k, i );
			if( mindisfrom[i] < minscore ) // muscle
			{
				im = i;
				minscore = mindisfrom[i];
			}
		}
		jm = nearest[im];
		if( jm < im ) 
		{
			j=jm; jm=im; im=j;
		}


		prevnode = hist[im];
		nmemim = nmemar[im];
//		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		intpt = topol[k][0];
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		prevnode = hist[jm];
		nmemjm = nmemar[jm];
//		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		intpt = topol[k][1];
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );


		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

		mindisfrom[im] = 999.9;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
#if 0
                		tmpdouble = eff[miniim][maxiim] =
				MIN( eff0, eff1 ) * sueff1 + ( eff0 + eff1 ) * sueff05; 
#else
                tmpdouble = eff[miniim][maxiim] =
				(clusterfuncpt[0])( eff0, eff1 );
#endif
				if( tmpdouble < mindisfrom[i]  )
				{
					mindisfrom[i] = tmpdouble;
					nearest[i] = im;
				}
				if( tmpdouble < mindisfrom[im]  )
				{
					mindisfrom[im] = tmpdouble;
					nearest[im] = i;
				}
				if( nearest[i] == jm )
				{
					nearest[i] = im;
				}
            }
        }
#if 0
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		strcpy( tree[im], treetmp );
#else
		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );
#endif

		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
//		free( (void *)eff[jm] ); eff[jm] = NULL;

#if 1 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
				}
				if( eff[miniim][maxiim] > mindisfrom[i] )
					setnearest_double_fullmtx( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#endif


#if 0
        fprintf( stdout, "\nvSTEP-%03d:\n", k+1 );
		fprintf( stdout, "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][0][i]+1 );
        fprintf( stdout, "\n" );
		fprintf( stdout, "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][1][i]+1 );
        fprintf( stdout, "\n" );
#endif
    }
	fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
	fclose( fp );
#if 0
	FreeCharMtx( tree );
#else
	free( tree[0] );
	free( tree );
#endif
	free( treetmp );
	free( nametmp );
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	free( mindisfrom );
	free( nearest );
}

void fixed_supg_double_treeout_constrained( int nseq, double **eff, int ***topol, double **len, char **name, int ngroup, int **groups )
{
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	double tmpdouble;
	double eff1, eff0;
	static double *tmptmplen = NULL;
	static int *hist = NULL;
	static Bchain *ac = NULL;
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti, *acptj;
	int *pt1, *pt2, *pt11, *pt22;
	static int *nmemar;
	int nmemim, nmemjm;
	double minscore;
	int *nearest = NULL; // by D.Mathog, a guess
	double *mindisfrom = NULL; // by D.Mathog, a guess
	static char **tree;
	static char *treetmp;
	static char *nametmp, *nameptr, *tmpptr;
	FILE *fp;
	double (*clusterfuncpt[1])(double,double);
	char namec;
	int *testtopol, **inconsistent;
	int **inconsistentpairlist;
	int ninconsistentpairs;
	int *warned;
	int allinconsistent;
	int firsttime;

	increaseintergroupdistancesfullmtx( eff, ngroup, groups, nseq );

	sueff1 = 1 - sueff_global;
	sueff05 = sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_double;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_double;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_double;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}





#if 0
	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateDoubleVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateDoubleVec( njob );
		nearest = AllocateIntVec( njob );
		treetmp = AllocateCharVec( njob*150 );
		nametmp = AllocateCharVec( 91 );
		tree = AllocateCharMtx( njob, njob*150 );
	}
	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<90; j++ ) nametmp[j] = 0;
		for( j=0; j<90; j++ ) 
		{
			if( name[i][j] == 0 )
				break;
			else if( isalnum( name[i][j] ) )
				nametmp[j] = name[i][j];
			else
				nametmp[j] = '_';
		}
		nametmp[90] = 0;
//		sprintf( tree[i], "%d_%.60s", i+1, nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		sprintf( tree[i], "\n%d_%.60s\n", i+1, nameptr );
	}

#else

	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateDoubleVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateDoubleVec( njob );
		nearest = AllocateIntVec( njob );
//		treetmp = AllocateCharVec( njob * ( B + 100 ) ); // nagasugi?
		treetmp = NULL; // kentou 2013/06/12
		nametmp = AllocateCharVec( 1000 ); // nagasugi
//		tree = AllocateCharMtx( njob, njob*600 );
		tree = AllocateCharMtx( njob, 0 );
		testtopol = AllocateIntVec( njob + 1 );
		inconsistent = AllocateIntMtx( njob, njob ); // muda
		inconsistentpairlist = AllocateIntMtx( njob*(njob-1)/2+1, 2 ); // muda
		warned = AllocateIntVec( ngroup );
	}

	
	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}

#endif








	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

	for( i=0; i<nseq; i++ ) setnearest_double_fullmtx( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	reporterr(       "\n" );
	ninconsistentpairs = 0;
	for( k=0; k<nseq-1; k++ )
	{
		if( k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );



//		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) for( acptj=acpti->next; acptj!=NULL; acptj=acptj->next ) inconsistent[acpti->pos][acptj->pos] = 0;
		for( i=0; i<ninconsistentpairs; i++ ) inconsistent[inconsistentpairlist[i][0]][inconsistentpairlist[i][1]] = 0;
		ninconsistentpairs = 0;
		firsttime = 1;
		while( 1 )
		{
			if( firsttime )
			{
				firsttime = 0;
				minscore = 999.9;
				for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
				{
					i = acpti->pos;
//					reporterr(       "k=%d i=%d\n", k, i );
					if( mindisfrom[i] < minscore ) // muscle
					{
						im = i;
						minscore = mindisfrom[i];
					}
				}
				jm = nearest[im];
				if( jm < im ) 
				{
					j=jm; jm=im; im=j;
				}
			}
			else
			{
				minscore = 999.9;
				for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
				{
					i = acpti->pos;
//					reporterr(       "k=%d i=%d\n", k, i );
					for( acptj=acpti->next; acptj!=NULL; acptj=acptj->next )
					{
						j = acptj->pos;
						if( !inconsistent[i][j] && (tmpdouble=eff[i][j]) < minscore )
						{
							minscore = tmpdouble;
							im = i; jm = j;
						}
					}
					for( acptj=ac; (acptj&&acptj->pos!=i); acptj=acptj->next )
					{
						j = acptj->pos;
						if( !inconsistent[j][i] && (tmpdouble=eff[j][i]) < minscore )
						{
							minscore = tmpdouble;
							im = j; jm = i;
						}
					}
				}
			}

			allinconsistent = 1;
			for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
			{
				for( acptj=acpti->next; acptj!=NULL; acptj=acptj->next )
				{
					if( inconsistent[acpti->pos][acptj->pos] == 0 )
					{
						allinconsistent = 0;
						goto exitloop_d;
					}
				}
			}
			exitloop_d:

			if( allinconsistent )
			{
				reporterr(       "\n\n\nPlease check whether the grouping is possible.\n\n\n" );
				exit( 1 );
			}
#if 1
			intpt = testtopol;
			prevnode = hist[im];
			if( prevnode == -1 )
			{
				*intpt++ = im;
			}
			else
			{
				for( intpt2=topol[prevnode][0]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
				for( intpt2=topol[prevnode][1]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
			}
			
			prevnode = hist[jm];
			if( prevnode == -1 )
			{
				*intpt++ = jm;
			}
			else
			{
				for( intpt2=topol[prevnode][0]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
				for( intpt2=topol[prevnode][1]; *intpt2!=-1; )
					*intpt++ = *intpt2++;
			}
			*intpt = -1;
//			reporterr(       "testtopol = \n" );
//       	for( i=0; testtopol[i]>-1; i++ ) reporterr(       " %03d", testtopol[i]+1 );
//			reporterr(       "\n" );
#endif
			for( i=0; i<ngroup; i++ )
			{
//				reporterr(       "groups[%d] = \n", i );
//  		     	for( j=0; groups[i][j]>-1; j++ ) reporterr(       " %03d", groups[i][j]+1 );
//				reporterr(       "\n" );
				if( overlapmember( testtopol, groups[i] ) )
				{
					if( !includemember( testtopol, groups[i] ) && !includemember( groups[i], testtopol ) )
					{
						if( !warned[i] )
						{
							warned[i] = 1;
							reporterr(       "\n###################################################################\n" );
							reporterr(       "# WARNING: Group %d is forced to be a monophyletic cluster.\n", i+1 );
							reporterr(       "###################################################################\n" );
						}
						inconsistent[im][jm] = 1;
						inconsistentpairlist[ninconsistentpairs][0] = im;
						inconsistentpairlist[ninconsistentpairs][1] = jm;
						ninconsistentpairs++;
						break;
					}
				}
			}
			if( i == ngroup )
			{
//				reporterr(       "OK\n" );
				break;
			}
		}






		prevnode = hist[im];
		nmemim = nmemar[im];
//		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		intpt = topol[k][0];
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		prevnode = hist[jm];
		nmemjm = nmemar[jm];
//		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		intpt = topol[k][1];
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );
		if( len[k][0] < 0.0 ) len[k][0] = 0.0;
		if( len[k][1] < 0.0 ) len[k][1] = 0.0;


		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

		mindisfrom[im] = 999.9;
		eff[im][jm] = 999.9;
//		eff[im][jm-im] = 999.9; // bug??

		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
#if 0
                		tmpdouble = eff[miniim][maxiim] =
				MIN( eff0, eff1 ) * sueff1 + ( eff0 + eff1 ) * sueff05; 
#else
                tmpdouble = eff[miniim][maxiim] =
				(clusterfuncpt[0])( eff0, eff1 );
#endif

#if 1
				if( tmpdouble < mindisfrom[i]  )
				{
					mindisfrom[i] = tmpdouble;
					nearest[i] = im;
				}
				if( tmpdouble < mindisfrom[im]  )
				{
					mindisfrom[im] = tmpdouble;
					nearest[im] = i;
				}
				if( nearest[i] == jm )
				{
					nearest[i] = im;
				}
#endif
            }
        }
#if 0
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		strcpy( tree[im], treetmp );
#else
		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );
#endif

		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
//		free( (void *)eff[jm] ); eff[jm] = NULL;

#if 1 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
				}
				if( eff[miniim][maxiim] > mindisfrom[i] )
					setnearest_double_fullmtx( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#endif


#if 0
        fprintf( stdout, "\ncSTEP-%03d:\n", k+1 );
		fprintf( stdout, "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][0][i]+1 );
        fprintf( stdout, "\n" );
		fprintf( stdout, "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][1][i]+1 );
        fprintf( stdout, "\n" );
#endif
    }
	fp = fopen( "infile.tree", "w" );
		fprintf( fp, "%s\n", treetmp );
	fclose( fp );
#if 0
	FreeCharMtx( tree );
#else
	free( tree[0] );
	free( tree );
#endif
	free( treetmp );
	free( nametmp );
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	free( mindisfrom );
	free( nearest );
	free( testtopol );
	FreeIntMtx( inconsistent );
	FreeIntMtx( inconsistentpairlist );
	free( warned );
}

void fixed_musclesupg_double_realloc_nobk_halfmtx_memsave( int nseq, double **eff, int ***topol, double **len, Treedep *dep, int progressout, int efffree )
{
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt;
	double tmpdouble;
	double eff1, eff0;
	double *tmptmplen = NULL; // static TLS -> local, 2012/02/25
	int *hist = NULL; // static TLS -> local, 2012/02/25 
	Bchain *ac = NULL; // static TLS -> local, 2012/02/25
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti;
	int *pt1, *pt2, *pt11;
	int *nmemar; // static TLS -> local, 2012/02/25
	int nmemim, nmemjm;
	double minscore;
	int *nearest = NULL; // by Mathog, a guess
	double *mindisfrom = NULL; // by Mathog, a guess
	double (*clusterfuncpt[1])(double,double);


	sueff1 = 1 - (double)sueff_global;
	sueff05 = (double)sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_double;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_double;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_double;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}

	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateFloatVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateFloatVec( njob );
		nearest = AllocateIntVec( njob );
	}

	
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

	for( i=0; i<nseq; i++ ) setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	if( progressout ) reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( progressout && k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );

		minscore = 999.9;
		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
		{
			i = acpti->pos;
//			reporterr(       "k=%d i=%d\n", k, i );
			if( mindisfrom[i] < minscore ) // muscle
			{
				im = i;
				minscore = mindisfrom[i];
			}
		}
		jm = nearest[im];
		if( jm < im ) 
		{
			j=jm; jm=im; im=j;
		}


		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];
		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( 2 ) * sizeof( int ) ); // memsave
//		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
//				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
//				pt22 = pt2;
			}
#if 1
			*intpt++ = *pt11;
			*intpt = -1;
#else
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
#endif
		}

		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;
		nmemjm = nmemar[jm];
		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( 2 ) * sizeof( int ) );
//		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
//				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
//				pt22 = pt2;
			}
#if 1
			*intpt++ = *pt11;
			*intpt = -1;
#else
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
#endif
		}

		minscore *= 0.5;

		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );

		if( dep ) dep[k].distfromtip = minscore;

		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

		mindisfrom[im] = 999.9;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim-miniim];
				eff1 = eff[minijm][maxijm-minijm];
				tmpdouble = eff[miniim][maxiim-miniim] =
#if 0
				MIN( eff0, eff1 ) * sueff1 + ( eff0 + eff1 ) * sueff05; 
#else
				(clusterfuncpt[0])( eff0, eff1 );
#endif
				if( tmpdouble < mindisfrom[i]  )
				{
					mindisfrom[i] = tmpdouble;
					nearest[i] = im;
				}
				if( tmpdouble < mindisfrom[im]  )
				{
					mindisfrom[im] = tmpdouble;
					nearest[im] = i;
				}
				if( nearest[i] == jm )
				{
					nearest[i] = im;
				}
            }
        }

//		reporterr(       "im,jm=%d,%d\n", im, jm );
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
		if( efffree )
		{
			free( (void *)eff[jm] ); eff[jm] = NULL;
		}

#if 1 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
				}
				if( eff[miniim][maxiim-miniim] > mindisfrom[i] )
					setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#endif


#if 0
        fprintf( stdout, "vSTEP-%03d:\n", k+1 );
		fprintf( stdout, "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][0][i]+1 );
        fprintf( stdout, "\n" );
		fprintf( stdout, "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][1][i]+1 );
        fprintf( stdout, "\n" );
#endif
    }
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	free( mindisfrom );
	free( nearest );
}
void fixed_musclesupg_double_realloc_nobk_halfmtx( int nseq, double **eff, int ***topol, double **len, Treedep *dep, int progressout, int efffree )
{
	int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	double tmpdouble;
	double eff1, eff0;
	double *tmptmplen = NULL; // static TLS -> local, 2012/02/25
	int *hist = NULL; // static TLS -> local, 2012/02/25 
	Bchain *ac = NULL; // static TLS -> local, 2012/02/25
	int im = -1, jm = -1;
	Bchain *acjmnext, *acjmprev;
	int prevnode;
	Bchain *acpti;
	int *pt1, *pt2, *pt11, *pt22;
	int *nmemar; // static TLS -> local, 2012/02/25
	int nmemim, nmemjm;
	double minscore;
	int *nearest = NULL; // by Mathog, a guess
	double *mindisfrom = NULL; // by Mathog, a guess
	double (*clusterfuncpt[1])(double,double);


	sueff1 = 1 - (double)sueff_global;
	sueff05 = (double)sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_double;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_double;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_double;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}

	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = AllocateFloatVec( njob );
		ac = (Bchain *)malloc( njob * sizeof( Bchain ) );
		nmemar = AllocateIntVec( njob );
		mindisfrom = AllocateFloatVec( njob );
		nearest = AllocateIntVec( njob );
	}

	
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[nseq-1].next = NULL;

	for( i=0; i<nseq; i++ ) setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i ); // muscle

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		hist[i] = -1;
		nmemar[i] = 1;
	}

	if( progressout ) reporterr(       "\n" );
	for( k=0; k<nseq-1; k++ )
	{
		if( progressout && k % 10 == 0 ) reporterr(       "\r% 5d / %d", k, nseq );

		minscore = 999.9;
		for( acpti=ac; acpti->next!=NULL; acpti=acpti->next ) 
		{
			i = acpti->pos;
//			reporterr(       "k=%d i=%d\n", k, i );
			if( mindisfrom[i] < minscore ) // muscle
			{
				im = i;
				minscore = mindisfrom[i];
			}
		}
		jm = nearest[im];
		if( jm < im ) 
		{
			j=jm; jm=im; im=j;
		}


		prevnode = hist[im];
		if( dep ) dep[k].child0 = prevnode;
		nmemim = nmemar[im];
		intpt = topol[k][0] = (int *)realloc( topol[k][0], ( nmemim + 1 ) * sizeof( int ) );
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		prevnode = hist[jm];
		if( dep ) dep[k].child1 = prevnode;
		nmemjm = nmemar[jm];
		intpt = topol[k][1] = (int *)realloc( topol[k][1], ( nmemjm + 1 ) * sizeof( int ) );
		if( !intpt )
		{
			reporterr(       "Cannot reallocate topol\n" );
			exit( 1 );
		}
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = ( minscore - tmptmplen[im] );
		len[k][1] = ( minscore - tmptmplen[jm] );

		if( dep ) dep[k].distfromtip = minscore;

		tmptmplen[im] = minscore;

		hist[im] = k;
		nmemar[im] = nmemim + nmemjm;

		mindisfrom[im] = 999.9;
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
        {
			i = acpti->pos;
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim-miniim];
				eff1 = eff[minijm][maxijm-minijm];
				tmpdouble = eff[miniim][maxiim-miniim] =
#if 0
				MIN( eff0, eff1 ) * sueff1 + ( eff0 + eff1 ) * sueff05; 
#else
				(clusterfuncpt[0])( eff0, eff1 );
#endif
				if( tmpdouble < mindisfrom[i]  )
				{
					mindisfrom[i] = tmpdouble;
					nearest[i] = im;
				}
				if( tmpdouble < mindisfrom[im]  )
				{
					mindisfrom[im] = tmpdouble;
					nearest[im] = i;
				}
				if( nearest[i] == jm )
				{
					nearest[i] = im;
				}
            }
        }

//		reporterr(       "im,jm=%d,%d\n", im, jm );
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		acjmprev->next = acjmnext;
		if( acjmnext != NULL )
			acjmnext->prev = acjmprev;
		if( efffree )
		{
			free( (void *)eff[jm] ); eff[jm] = NULL;
		}

#if 1 // muscle seems to miss this.
		for( acpti=ac; acpti!=NULL; acpti=acpti->next )
		{
			i = acpti->pos;
			if( nearest[i] == im ) 
			{
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
				}
				if( eff[miniim][maxiim-miniim] > mindisfrom[i] )
					setnearest( nseq, ac, eff, mindisfrom+i, nearest+i, i );
			}
		}
#endif


#if 0
        fprintf( stdout, "vSTEP-%03d:\n", k+1 );
		fprintf( stdout, "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][0][i]+1 );
        fprintf( stdout, "\n" );
		fprintf( stdout, "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][1][i]+1 );
        fprintf( stdout, "\n" );
#endif
    }
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	free( (void *)nmemar ); nmemar = NULL;
	free( mindisfrom );
	free( nearest );
}









void veryfastsupg_double_loadtree( int nseq, double **eff, int ***topol, double **len, char **name )
{
    int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	double eff1, eff0;
    int *hist = NULL;
	Achain *ac = NULL;
	double minscore;
	char **tree;
	char *treetmp;
	int im = -1, jm = -1;
	int prevnode, acjmnext, acjmprev;
	int *pt1, *pt2, *pt11, *pt22;
	FILE *fp;
	int node[2];
	double lenfl[2];
	char *nametmp, *nameptr, *tmpptr; //static?
	char namec;

	fp = fopen( "_guidetree", "r" );
	if( !fp )
	{
		reporterr(       "cannot open _guidetree\n" );
		exit( 1 );
	}


	if( !hist )
	{
//		treetmp = AllocateCharVec( njob*50 );
		treetmp = NULL;
//		tree = AllocateCharMtx( njob, njob*50 );
		tree = AllocateCharMtx( njob, 0 );
		nametmp = AllocateCharVec( 1000 ); // nagasugi
		hist = AllocateIntVec( njob );
		ac = (Achain *)malloc( njob * sizeof( Achain ) );
	}

	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<999; j++ ) nametmp[j] = 0;
		for( j=0; j<999; j++ ) 
		{
			namec = name[i][j];
			if( namec == 0 )
				break;
			else if( isalnum( namec ) || namec == '/' || namec == '=' || namec == '-' || namec == '{' || namec == '}' )
				nametmp[j] = namec;
			else
				nametmp[j] = '_';
		}
		nametmp[j] = 0;
//		sprintf( tree[i], "%d_l=%d_%.20s", i+1, nlen[i], nametmp+1 );
		if( outnumber )
			nameptr = strstr( nametmp, "_numo_e" ) + 8;
		else
			nameptr = nametmp + 1;

		if( (tmpptr=strstr( nameptr, "_oe_" )) ) nameptr = tmpptr + 4; // = -> _ no tame

		tree[i] = calloc( strlen( nametmp )+100, sizeof( char ) ); // suuji no bun de +100
		if( tree[i] == NULL )
		{
			reporterr(       "Cannot allocate tree!\n" );
			exit( 1 );
		}
		sprintf( tree[i], "\n%d_%.900s\n", i+1, nameptr );
	}
	
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = i+1;
		ac[i].prev = i-1;
//		ac[i].curr = i;
	}
	ac[nseq-1].next = -1;

    for( i=0; i<nseq; i++ ) hist[i] = -1;

	reporterr(       "\n" );
    for( k=0; k<nseq-1; k++ )
    {
		if( k % 10 == 0 ) reporterr(       "%d / %d\r", k, nseq );

#if 0
		minscore = 99999.9;
		for( i=0; ac[i].next!=-1; i=ac[i].next ) 
		{
			for( j=ac[i].next; j!=-1; j=ac[j].next )
	        {
				tmpdouble = eff[i][j];
				if( tmpdouble < minscore )
				{
					minscore = tmpdouble;
					im = i; jm = j;
				}
			}
		}
#else
		lenfl[0] = lenfl[1] = -1.0;
		loadtreeoneline( node, lenfl, fp );
		im = node[0];
		jm = node[1];
		minscore = eff[im][jm];

		if( im > nseq-1 || jm > nseq-1 || tree[im] == NULL || tree[jm] == NULL )
		{
			reporterr(       "\n\nCheck the guide tree.\n" );
			reporterr(       "im=%d, jm=%d\n", im+1, jm+1 );
			reporterr(       "Please use newick2mafft.rb to generate a tree file from a newick tree.\n\n" );
			exit( 1 );
		}


//		reporterr(       "im=%d, jm=%d, minscore = %f\n", im, jm, minscore );


		if( lenfl[0] == -1.0 || lenfl[1] == -1.0 )
		{
			reporterr(       "\n\nWARNING: Branch length is not given.\n" );
			exit( 1 );
		}

		if( lenfl[0] < 0.0 ) lenfl[0] = 0.0;
		if( lenfl[1] < 0.0 ) lenfl[1] = 0.0;
#endif

//		reporterr(       "im=%d, jm=%d\n", im, jm );

		intpt = topol[k][0];
		prevnode = hist[im];
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		intpt = topol[k][1];
		prevnode = hist[jm];
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

#if 0
		len[k][0] = minscore - tmptmplen[im];
		len[k][1] = minscore - tmptmplen[jm];
#else
		len[k][0] = lenfl[0];
		len[k][1] = lenfl[1];
#endif


		hist[im] = k;

		for( i=0; i!=-1; i=ac[i].next )
        {
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                eff[miniim][maxiim] =
                MIN( eff0, eff1 ) * ( 1.0 - sueff_global ) +
				( eff0 + eff1 ) * 0.5 * sueff_global;
            }
        }
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		ac[acjmprev].next = acjmnext;
		if( acjmnext != -1 )
			ac[acjmnext].prev = acjmprev;


		treetmp = realloc( treetmp, strlen( tree[im] ) + strlen( tree[jm] ) + 100 ); // 22 de juubunn (:%7,:%7) %7 ha minus kamo
		if( !treetmp )
		{
			reporterr(       "Cannot allocate treetmp\n" );
			exit( 1 );
		}
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		free( tree[im] );
		free( tree[jm] );
		tree[im] = calloc( strlen( treetmp )+1, sizeof( char ) );
		tree[jm] = NULL;
		if( tree[im] == NULL )
		{
			reporterr(       "Cannot reallocate tree!\n" );
			exit( 1 );
		}
		strcpy( tree[im], treetmp );

//		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
//		strcpy( tree[im], treetmp );

#if 0
        fprintf( stdout, "STEP-%03d:\n", k+1 );
		fprintf( stdout, "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][0][i] );
        fprintf( stdout, "\n" );
		fprintf( stdout, "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][1][i] );
        fprintf( stdout, "\n" );
#endif
    }
	fclose( fp );


	fp = fopen( "infile.tree", "w" );
	fprintf( fp, "%s\n", treetmp );
//	fprintf( fp, "by veryfastsupg_double_loadtree\n" );
	fclose( fp );

#if 1
	reporterr(       "\n" );
	free( hist );
	free( (char *)ac );
	FreeCharMtx( tree );
	free( treetmp );
	free( nametmp );
#endif

#if 0
//	reporterr(       "reconstructing eff[][]\n" ); // Tsune ni hat2 ha aru node koreha iranai.
    for( k=0; k<nseq; k++ ) for( i=0; i<nseq; i++ ) eff[i][k] = 0.0;
    for( k=0; k<nseq-1; k++ )
	{
		reporterr(       "len[k][0], len[k][1] = %f, %f\n", len[k][0], len[k][1] );
        for( i=0; (im=topol[k][0][i])>-1; i++ )
		{
			reporterr(       " %03d", im );
		}
		fprintf( stdout, "\n" );
        for( i=0; (jm=topol[k][1][i])>-1; i++ ) 
		{
			reporterr(       " %03d", jm );
		}
        for( i=0; (im=topol[k][0][i])>-1; i++ ) for( j=0; (jm=topol[k][1][j])>-1; j++ ) 
		{
			eff[im][jm] += len[k][0] + len[k][1];
			eff[jm][im] += len[k][0] + len[k][1];
		}
	}
#endif
}

#if 0
void veryfastsupg_double( int nseq, double **eff, int ***topol, double **len )
{
    int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	double tmpdouble;
	double eff1, eff0;
	static double *tmptmplen = NULL;
    static int *hist = NULL;
	static Achain *ac = NULL;
	double minscore;
	int im = -1, jm = -1;
	int prevnode, acjmnext, acjmprev;
	int *pt1, *pt2, *pt11, *pt22;
	if( !hist )
	{
		hist = AllocateIntVec( njob );
		tmptmplen = (double *)malloc( njob * sizeof( double ) );
		ac = (Achain *)malloc( njob * sizeof( Achain ) );
	}
	
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = i+1;
		ac[i].prev = i-1;
//		ac[i].curr = i;
	}
	ac[nseq-1].next = -1;

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) hist[i] = -1;

	reporterr(       "\n" );
    for( k=0; k<nseq-1; k++ )
    {
		if( k % 10 == 0 ) reporterr(       "%d / %d\r", k, nseq );

		minscore = 99999.9;
		for( i=0; ac[i].next!=-1; i=ac[i].next ) 
		{
			for( j=ac[i].next; j!=-1; j=ac[j].next )
	        {
				tmpdouble = eff[i][j];
				if( tmpdouble < minscore )
				{
					minscore = tmpdouble;
					im = i; jm = j;
				}
			}
		}

//		reporterr(       "im=%d, jm=%d\n", im, jm );

		intpt = topol[k][0];
		prevnode = hist[im];
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		intpt = topol[k][1];
		prevnode = hist[jm];
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = minscore - tmptmplen[im];
		len[k][1] = minscore - tmptmplen[jm];

		tmptmplen[im] = minscore;

		hist[im] = k;

		for( i=0; i!=-1; i=ac[i].next )
        {
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                eff[miniim][maxiim] =
                MIN( eff0, eff1 ) * ( 1.0 - sueff_global ) +
				( eff0 + eff1 ) * 0.5 * sueff_global;
            }
        }
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		ac[acjmprev].next = acjmnext;
		if( acjmnext != -1 )
			ac[acjmnext].prev = acjmprev;
#if 0
        fprintf( stdout, "STEP-%03d:\n", k+1 );
		fprintf( stdout, "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][0][i] );
        fprintf( stdout, "\n" );
		fprintf( stdout, "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][1][i] );
        fprintf( stdout, "\n" );
#endif
    }
#if 1
	reporterr(       "\n" );
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
#endif
}
#endif

void veryfastsupg_double_outtree( int nseq, double **eff, int ***topol, double **len, char **name ) // not used
{
    int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	double tmpdouble;
	double eff1, eff0;
	static double *tmptmplen = NULL;
    static int *hist = NULL;
	static Achain *ac = NULL;
	double minscore;
	static char **tree;
	static char *treetmp;
	static char *nametmp;
	FILE *fpout;
	int im = -1, jm = -1;
	int prevnode, acjmnext, acjmprev;
	int *pt1, *pt2, *pt11, *pt22;
	double (*clusterfuncpt[1])(double,double);


	sueff1 = 1 - sueff_global;
	sueff05 = sueff_global * 0.5;
	if ( treemethod == 'X' )
		clusterfuncpt[0] = cluster_mix_double;
	else if ( treemethod == 'E' )
		clusterfuncpt[0] = cluster_average_double;
	else if ( treemethod == 'q' )
		clusterfuncpt[0] = cluster_minimum_double;
	else
	{
		reporterr(       "Unknown treemethod, %c\n", treemethod );
		exit( 1 );
	}

	if( !hist )
	{
		treetmp = AllocateCharVec( njob*50 );
		tree = AllocateCharMtx( njob, njob*50 );
		hist = AllocateIntVec( njob );
		tmptmplen = (double *)malloc( njob * sizeof( double ) );
		ac = (Achain *)malloc( njob * sizeof( Achain ) );
		nametmp = AllocateCharVec( 31 );
	}

//	for( i=0; i<nseq; i++ ) sprintf( tree[i], "%d", i+1 );
    for( i=0; i<nseq; i++ )
	{
		for( j=0; j<30; j++ ) nametmp[j] = 0;
		for( j=0; j<30; j++ ) 
		{
			if( isalnum( name[i][j] ) )
				nametmp[j] = name[i][j];
			else
				nametmp[j] = '_';
		}
		nametmp[30] = 0;
		sprintf( tree[i], "%d_%.20s", i+1, nametmp+1 );
	}
	
	for( i=0; i<nseq; i++ )
	{
		ac[i].next = i+1;
		ac[i].prev = i-1;
//		ac[i].curr = i;
	}
	ac[nseq-1].next = -1;

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) hist[i] = -1;

	reporterr(       "\n" );
    for( k=0; k<nseq-1; k++ )
    {
		if( k % 10 == 0 ) reporterr(       "%d / %d\r", k, nseq );

		minscore = 99999.9;
		for( i=0; ac[i].next!=-1; i=ac[i].next ) 
		{
			for( j=ac[i].next; j!=-1; j=ac[j].next )
	        {
				tmpdouble = eff[i][j];
				if( tmpdouble < minscore )
				{
					minscore = tmpdouble;
					im = i; jm = j;
				}
			}
		}

//		reporterr(       "im=%d, jm=%d\n", im, jm );

		intpt = topol[k][0];
		prevnode = hist[im];
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		intpt = topol[k][1];
		prevnode = hist[jm];
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		minscore *= 0.5;

		len[k][0] = minscore - tmptmplen[im];
		len[k][1] = minscore - tmptmplen[jm];

		tmptmplen[im] = minscore;

		hist[im] = k;

		for( i=0; i!=-1; i=ac[i].next )
        {
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                eff[miniim][maxiim] =
				(clusterfuncpt[0])( eff0, eff1 );
            }
        }
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		ac[acjmprev].next = acjmnext;
		if( acjmnext != -1 )
			ac[acjmnext].prev = acjmprev;

		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], len[k][0], tree[jm], len[k][1] );
		strcpy( tree[im], treetmp );
#if 0
        fprintf( stdout, "STEP-%03d:\n", k+1 );
		fprintf( stdout, "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][0][i] );
        fprintf( stdout, "\n" );
		fprintf( stdout, "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][1][i] );
        fprintf( stdout, "\n" );
#endif
    }
	fpout = fopen( "infile.tree", "w" );
	fprintf( fpout, "%s\n", treetmp );
//	fprintf( fpout, "by veryfastsupg_double_outtree\n" );
	fclose( fpout );
#if 1
	reporterr(       "\n" );
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
	FreeCharMtx( tree );
	free( treetmp );
	free( nametmp );
#endif
}

void veryfastsupg( int nseq, double **oeff, int ***topol, double **len )
{
    int i, j, k, miniim, maxiim, minijm, maxijm;
	int *intpt, *intpt2;
	int tmpint;
	int eff1, eff0;
	static double *tmptmplen = NULL;
	static int **eff = NULL;
    static int *hist = NULL;
	static Achain *ac = NULL;
	int minscore;
	double minscoref;
	int im = -1, jm = -1;
	int prevnode, acjmnext, acjmprev;
	int *pt1, *pt2, *pt11, *pt22;
	if( !eff )
	{
		eff = AllocateIntMtx( njob, njob );
		hist = AllocateIntVec( njob );
		tmptmplen = (double *)malloc( njob * sizeof( double ) );
		ac = (Achain *)malloc( njob * sizeof( Achain ) );
	}
	
	for( i=0; i<nseq; i++ ) 
	{
		for( j=0; j<nseq; j++ ) 
		{
			eff[i][j] = (int)( oeff[i][j] * INTMTXSCALE + 0.5 );
		}
	}

	for( i=0; i<nseq; i++ )
	{
		ac[i].next = i+1;
		ac[i].prev = i-1;
//		ac[i].curr = i;
	}
	ac[nseq-1].next = -1;

	for( i=0; i<nseq; i++ ) tmptmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) hist[i] = -1;

	reporterr(       "\n" );
    for( k=0; k<nseq-1; k++ )
    {
		if( k % 10 == 0 ) reporterr(       "%d / %d\r", k, nseq );

		minscore = INTMTXSCALE*4;
		for( i=0; ac[i].next!=-1; i=ac[i].next ) 
		{
			for( j=ac[i].next; j!=-1; j=ac[j].next )
	        {
				tmpint = eff[i][j];
				if( tmpint < minscore )
				{
					minscore = tmpint;
					im = i; jm = j;
				}
			}
		}
		minscoref = (double)minscore * 0.5 / ( INTMTXSCALE );

//		reporterr(       "im=%d, jm=%d\n", im, jm );

#if 1
		intpt = topol[k][0];
		prevnode = hist[im];
		if( prevnode == -1 )
		{
			*intpt++ = im;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}

		intpt = topol[k][1];
		prevnode = hist[jm];
		if( prevnode == -1 )
		{
			*intpt++ = jm;
			*intpt = -1;
		}
		else
		{
			pt1 = topol[prevnode][0];
			pt2 = topol[prevnode][1];
			if( *pt1 > *pt2 )
			{
				pt11 = pt2;
				pt22 = pt1;
			}
			else
			{
				pt11 = pt1;
				pt22 = pt2;
			}
			for( intpt2=pt11; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			for( intpt2=pt22; *intpt2!=-1; )
				*intpt++ = *intpt2++;
			*intpt = -1;
		}
#else
		intpt = topol[k][0];
        for( i=0; i<nseq; i++ )
            if( pair[im][i] > -2 )
				*intpt++ = i;
		*intpt = -1;

		intpt = topol[k][1];
        for( i=0; i<nseq; i++ )
            if( pair[jm][i] > -2 )
				*intpt++ = i;
		*intpt = -1;
#endif

		len[k][0] = minscoref - tmptmplen[im];
		len[k][1] = minscoref - tmptmplen[jm];

		tmptmplen[im] = minscoref;

		hist[im] = k;

		for( i=0; i!=-1; i=ac[i].next )
        {
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                eff[miniim][maxiim] =
                MIN( eff0, eff1 ) * ( 1.0 - sueff_global ) + // int??
				( eff0 + eff1 ) * 0.5 * sueff_global;        // int??
            }
        }
		acjmprev = ac[jm].prev; 
		acjmnext = ac[jm].next; 
		ac[acjmprev].next = acjmnext;
		if( acjmnext != -1 )
			ac[acjmnext].prev = acjmprev;
#if 0
        fprintf( stdout, "STEP-%03d:\n", k+1 );
		fprintf( stdout, "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][0][i] );
        fprintf( stdout, "\n" );
		fprintf( stdout, "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) fprintf( stdout, " %03d", topol[k][1][i] );
        fprintf( stdout, "\n" );
#endif
    }
#if 1
	FreeIntMtx( eff ); eff = NULL;
	free( (void *)tmptmplen ); tmptmplen = NULL;
	free( hist ); hist = NULL;
	free( (char *)ac ); ac = NULL;
#endif
}

void fastsupg( int nseq, double **oeff, int ***topol, double **len )
{
    int i, j, k, miniim, maxiim, minijm, maxijm;
#if 0
	double eff[nseq][nseq];
    char pair[njob][njob];
#else
	static double *tmplen;
	int *intpt;
	double tmpdouble;
	double eff1, eff0;
	static double **eff = NULL;
    static char **pair = NULL;
	static Achain *ac;
	double minscore;
	int im = -1, jm = -1;
	if( !eff )
	{
		eff = AllocateFloatMtx( njob, njob );
		pair = AllocateCharMtx( njob, njob );
		tmplen = AllocateFloatVec( njob );
		ac = (Achain *)calloc( njob, sizeof( Achain ) );
	}
#endif
	
	for( i=0; i<nseq; i++ ) 
	{
		for( j=0; j<nseq; j++ ) 
		{
			eff[i][j] = (double)oeff[i][j];
		}
	}

	for( i=0; i<nseq; i++ )
	{
		ac[i].next = i+1;
		ac[i].prev = i-1;
//		ac[i].curr = i;
	}
	ac[nseq-1].next = -1;

	for( i=0; i<nseq; i++ ) tmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) pair[i][j] = 0;
    for( i=0; i<nseq; i++ ) pair[i][i] = 1;

	reporterr(       "\n" );
    for( k=0; k<nseq-1; k++ )
    {
		if( k % 10 == 0 ) reporterr(       "%d / %d\r", k, nseq );

		minscore = 9999.0;
		for( i=0; ac[i].next!=-1; i=ac[i].next ) 
//		for( i=0; i<nseq-1; i++ ) 
		{
			for( j=ac[i].next; j!=-1; j=ac[j].next )
//			for( j=i+1; j<nseq; j++ ) 
	        {
				tmpdouble = eff[i][j];
				if( tmpdouble < minscore )
				{
					minscore = tmpdouble;
					im = i; jm = j;
				}
			}
		}

//		reporterr(       "im=%d, jm=%d\n", im, jm );

		intpt = topol[k][0];
        for( i=0; i<nseq; i++ )
            if( pair[im][i] > 0 )
				*intpt++ = i;
		*intpt = -1;

		intpt = topol[k][1];
        for( i=0; i<nseq; i++ )
            if( pair[jm][i] > 0 )
				*intpt++ = i;
		*intpt = -1;

		minscore /= 2.0;

		len[k][0] = (double)minscore - tmplen[im];
		len[k][1] = (double)minscore - tmplen[jm];

		tmplen[im] = (double)minscore;

        for( i=0; i<nseq; i++ ) pair[im][i] += ( pair[jm][i] > 0 );
        for( i=0; i<nseq; i++ ) pair[jm][i] = 0;

//		for( i=0; i<nseq; i++ )
		for( i=0; i!=-1; i=ac[i].next )
        {
            if( i != im && i != jm )
            {
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                eff[miniim][maxiim] =
                MIN( eff0, eff1 ) * ( 1.0 - sueff_global ) +
				( eff0 + eff1 ) * 0.5 * sueff_global;
//        		eff[minijm][maxijm] = 9999.0;
            }
        }
		ac[ac[jm].prev].next = ac[jm].next;
		ac[ac[jm].next].prev = ac[jm].prev;
//		eff[im][jm] = 9999.0;
#if 0
        reporterr(       "STEP-%03d:\n", k+1 );
		reporterr(       "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) reporterr(       " %03d", topol[k][0][i] );
        reporterr(       "\n" );
		reporterr(       "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) reporterr(       " %03d", topol[k][1][i] );
        reporterr(       "\n" );
#endif
    }
	reporterr(       "\n" );

//	FreeFloatMtx( eff );
//	FreeCharMtx( pair );
//	FreeFloatVec( tmplen );
//	free( ac );
}
void supg( int nseq, double **oeff, int ***topol, double **len )
{
    int i, j, k, miniim, maxiim, minijm, maxijm;
#if 0
	double eff[nseq][nseq];
    char pair[njob][njob];
#else
	static double *tmplen;
	int *intpt;
	double **doubleptpt;
	double *doublept;
	double tmpdouble;
	double eff1, eff0;
	static double **eff = NULL;
    static char **pair = NULL;
	if( !eff )
	{
		eff = AllocateFloatMtx( njob, njob );
		pair = AllocateCharMtx( njob, njob );
		tmplen = AllocateFloatVec( njob );
	}
#endif

	
	for( i=0; i<nseq; i++ ) 
	{
		for( j=0; j<nseq; j++ ) 
		{
			eff[i][j] = (double)oeff[i][j];
		}
	}
	for( i=0; i<nseq; i++ ) tmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) pair[i][j] = 0;
    for( i=0; i<nseq; i++ ) pair[i][i] = 1;

    for( k=0; k<nseq-1; k++ )
    {
        double minscore = 9999.0;
        int im = -1, jm = -1;


		doubleptpt = eff;
        for( i=0; i<nseq-1; i++ ) 
		{
			doublept = *doubleptpt++ + i + 1;
			for( j=i+1; j<nseq; j++ )
	        {
				tmpdouble = *doublept++;
				if( tmpdouble < minscore )
				{
					minscore = tmpdouble;
					im = i; jm = j;
				}
			}
		}
		intpt = topol[k][0];
        for( i=0; i<nseq; i++ )
            if( pair[im][i] > 0 )
				*intpt++ = i;
		*intpt = -1;

		intpt = topol[k][1];
        for( i=0; i<nseq; i++ )
            if( pair[jm][i] > 0 )
				*intpt++ = i;
		*intpt = -1;

		len[k][0] = (double)minscore / 2.0 - tmplen[im];
		len[k][1] = (double)minscore / 2.0 - tmplen[jm];

		tmplen[im] = (double)minscore / 2.0;

        for( i=0; i<nseq; i++ ) pair[im][i] += ( pair[jm][i] > 0 );
        for( i=0; i<nseq; i++ ) pair[jm][i] = 0;

        for( i=0; i<nseq; i++ )
        {
            if( i != im && i != jm )
            {
#if 1
				if( i < im )
				{
					 miniim = i;
					 maxiim = im;
					 minijm = i;
					 maxijm = jm;
				}
				else if( i < jm )
				{
					 miniim = im;
					 maxiim = i;
					 minijm = i;
					 maxijm = jm;
				}
				else
				{
					 miniim = im;
					 maxiim = i;
					 minijm = jm;
					 maxijm = i;
				}
#else
				miniim = MIN( i, im );
				maxiim = MAX( i, im );
				minijm = MIN( i, jm );
				maxijm = MAX( i, jm );
#endif
#if 1
				eff0 = eff[miniim][maxiim];
				eff1 = eff[minijm][maxijm];
                eff[miniim][maxiim] =
                MIN( eff0, eff1 ) * ( 1.0 - sueff_global ) +
				( eff0 + eff1 ) * 0.5 * sueff_global;
#else
                MIN( eff[miniim][maxiim], eff[minijm][maxijm] ) * ( 1.0 - sueff_global ) +
				( eff[miniim][maxiim] + eff[minijm][maxijm] ) * 0.5 * sueff_global;
#endif
                eff[minijm][maxijm] = 9999.0;
            	eff[im][jm] = 9999.0;
            }
        }
#if DEBUG
        printf( "STEP-%03d:\n", k+1 );
		printf( "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) printf( " %03d", topol[k][0][i] );
        printf( "\n" );
		printf( "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) printf( " %03d", topol[k][1][i] );
        printf( "\n" );
#endif
    }
}

void spg( int nseq, double **oeff, int ***topol, double **len )
{
    int i, j, k;
	double tmplen[M];
#if 0
	double eff[nseq][nseq];
    char pair[njob][njob];
#else
	double **eff = NULL;
    char **pair = NULL;
	if( !eff )
	{
		eff = AllocateDoubleMtx( njob, njob );
		pair = AllocateCharMtx( njob, njob );
	}
#endif
	
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) eff[i][j] = oeff[i][j];
	for( i=0; i<nseq; i++ ) tmplen[i] = 0.0;
    for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) pair[i][j] = 0;
    for( i=0; i<nseq; i++ ) pair[i][i] = 1;

    for( k=0; k<nseq-1; k++ )
    {
        double minscore = 9999.0;
        int im = -1, jm = -1;
        int count;

        for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
        {
            if( eff[i][j] < minscore )
            {
                minscore = eff[i][j];
                im = i; jm = j;
            }
        }
        for( i=0, count=0; i<nseq; i++ )
            if( pair[im][i] > 0 )
            {
                topol[k][0][count] = i;
                count++;
            }
        topol[k][0][count] = -1;
        for( i=0, count=0; i<nseq; i++ )
            if( pair[jm][i] > 0 )
            {
                topol[k][1][count] = i;
                count++;
            }
        topol[k][1][count] = -1;

		len[k][0] = minscore / 2.0 - tmplen[im];
		len[k][1] = minscore / 2.0 - tmplen[jm];

		tmplen[im] = minscore / 2.0;

        for( i=0; i<nseq; i++ ) pair[im][i] += ( pair[jm][i] > 0 );
        for( i=0; i<nseq; i++ ) pair[jm][i] = 0;

        for( i=0; i<nseq; i++ )
        {
            if( i != im && i != jm )
            {
                eff[MIN(i,im)][MAX(i,im)] =
                MIN( eff[MIN(i,im)][MAX(i,im)], eff[MIN(i,jm)][MAX(i,jm)] );
                eff[MIN(i,jm)][MAX(i,jm)] = 9999.0;
            }
            eff[im][jm] = 9999.0;
        }
#if DEBUG
        printf( "STEP-%03d:\n", k+1 );
		printf( "len0 = %f\n", len[k][0] );
        for( i=0; topol[k][0][i]>-1; i++ ) printf( " %03d", topol[k][0][i] );
        printf( "\n" );
		printf( "len1 = %f\n", len[k][1] );
        for( i=0; topol[k][1][i]>-1; i++ ) printf( " %03d", topol[k][1][i] );
        printf( "\n" );
#endif
    }
}

double ipower( double x, int n )    /* n > 0  */
{
    double r;

    r = 1;
    while( n != 0 )
    {
        if( n & 1 ) r *= x;
        x *= x; n >>= 1;
    }
    return( r );
}

void countnode( int nseq, int ***topol, double **node ) /* node[j][i] != node[i][j] */
{
    int i, j, k, s1, s2;
    static double rootnode[M];

    if( nseq-2 < 0 )
	{
		reporterr(       "Too few sequence for countnode: nseq = %d\n", nseq );
		exit( 1 );
    }

    for( i=0; i<nseq; i++ ) rootnode[i] = 0;
    for( i=0; i<nseq-2; i++ )
    {
        for( j=0; topol[i][0][j]>-1; j++ )
            rootnode[topol[i][0][j]]++;
        for( j=0; topol[i][1][j]>-1; j++ )
            rootnode[topol[i][1][j]]++;
        for( j=0; topol[i][0][j]>-1; j++ )
        {
            s1 = topol[i][0][j];
            for( k=0; topol[i][1][k]>-1; k++ )
            {
                s2 = topol[i][1][k];
                node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2] - 1;
            }
        }
    }
    for( j=0; topol[nseq-2][0][j]>-1; j++ )
    {
        s1 = topol[nseq-2][0][j];
        for( k=0; topol[nseq-2][1][k]>-1; k++ )
        {
            s2 = topol[nseq-2][1][k];
            node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2];
        }
    }
}

void countnode_int( int nseq, int ***topol, int **node )  /* node[i][j] == node[j][i] */
{
    int i, j, k, s1, s2;
    int rootnode[M];

    for( i=0; i<nseq; i++ ) rootnode[i] = 0;
    for( i=0; i<nseq-2; i++ )
    {
        for( j=0; topol[i][0][j]>-1; j++ )
            rootnode[topol[i][0][j]]++;
        for( j=0; topol[i][1][j]>-1; j++ )
            rootnode[topol[i][1][j]]++;
        for( j=0; topol[i][0][j]>-1; j++ )
        {
            s1 = topol[i][0][j];
            for( k=0; topol[i][1][k]>-1; k++ )
            {
                s2 = topol[i][1][k];
                node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2] - 1;
            }
        }
    }
    for( j=0; topol[nseq-2][0][j]>-1; j++ )
    {
        s1 = topol[nseq-2][0][j];
        for( k=0; topol[nseq-2][1][k]>-1; k++ )
        {
            s2 = topol[nseq-2][1][k];
            node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2];
        }
    }
	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ ) 
		node[j][i] = node[i][j];
#if DEBUG
	reporterr(       "node[][] in countnode_int" );
	for( i=0; i<nseq; i++ ) 
	{
		for( j=0; j<nseq; j++ ) 
		{
			reporterr(       "%#3d", node[i][j] );
		}
		reporterr(       "\n" );
	}
#endif
}

void counteff_simple_double( int nseq, int ***topol, double **len, double *node )
{
    int i, j, s1, s2;
	double total;
	static double rootnode[M];
	static double eff[M];

#if DEBUG
	for( i=0; i<nseq; i++ ){
		reporterr(       "len0 = %f\n", len[i][0] );
		reporterr(       "len1 = %f\n", len[i][1] );
	}
#endif
    for( i=0; i<nseq; i++ )
	{
		rootnode[i] = 0.0;
		eff[i] = 1.0;
/*
		rootnode[i] = 1.0;
*/
	}
   	for( i=0; i<nseq-1; i++ )
   	{
       	for( j=0; (s1=topol[i][0][j]) > -1; j++ )
		{
           	rootnode[s1] += (double)len[i][0] * eff[s1];
			eff[s1] *= 0.5;
/*
           	rootnode[s1] *= 0.5;
*/
			
		}
       	for( j=0; (s2=topol[i][1][j]) > -1; j++ )
		{
           	rootnode[s2] +=  (double)len[i][1] * eff[s2];
			eff[s2] *= 0.5;
/*
           	rootnode[s2] *= 0.5;
*/
				
		}
	}
	for( i=0; i<nseq; i++ ) 
	{
#if 1 /* 97.9.29 */
		rootnode[i] += GETA3;
#endif
#if 0
		reporterr(       "### rootnode for %d = %f\n", i, rootnode[i] );
#endif
	}
#if 1
	total = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		total += rootnode[i];
	}
#else
	total = 1.0;
#endif
		
	for( i=0; i<nseq; i++ ) 
	{
		node[i] = rootnode[i] / total;
	}

#if 0
	reporterr(       "weight array in counteff_simple\n" );
	for( i=0; i<nseq; i++ )
		reporterr(       "%f\n", node[i] );
	printf( "\n" );
	exit( 1 );
#endif
}

void counteff_simple_double_nostatic_memsave( int nseq, int ***topol, double **len, Treedep *dep, double *node )
{
    int i, j, s1, s2;
	double total;
	double *rootnode;
	double *eff;
	int **localmem;
	int posinmem;

	rootnode = AllocateDoubleVec( nseq );
	eff = AllocateDoubleVec( nseq );
	localmem = AllocateIntMtx( 2, nseq+1 );

	for( i=0; i<nseq; i++ ) // 2014/06/07, fu no eff wo sakeru.
	{
		if( len[i][0] < 0.0 ) 
		{
			reporterr( "WARNING: negative branch length %f, step %d-0\n", len[i][0], i );
			len[i][0] = 0.0;
		}
		if( len[i][1] < 0.0 ) 
		{
			reporterr( "WARNING: negative branch length %f, step %d-1\n", len[i][1], i );
			len[i][1] = 0.0;
		}
	}
#if DEBUG
	for( i=0; i<nseq-1; i++ )
	{
		reporterr( "\nstep %d, group 0\n", i );
		for( j=0; topol[i][0][j]!=-1; j++) reporterr( "%3d ", topol[i][0][j] );
		reporterr( "\n", i );
		reporterr( "step %d, group 1\n", i );
		for( j=0; topol[i][1][j]!=-1; j++) reporterr( "%3d ", topol[i][1][j] );
		reporterr( "\n", i );
		reporterr(       "len0 = %f\n", len[i][0] );
		reporterr(       "len1 = %f\n", len[i][1] );
	}
#endif
    for( i=0; i<nseq; i++ )
	{
		rootnode[i] = 0.0;
		eff[i] = 1.0;
/*
		rootnode[i] = 1.0;
*/
	}
   	for( i=0; i<nseq-1; i++ )
   	{
		localmem[0][0] = -1;
		posinmem = 0;
		topolorder( njob, localmem[0], &posinmem, topol, dep, i, 0 );
		localmem[1][0] = -1;
		posinmem = 0;
		topolorder( njob, localmem[1], &posinmem, topol, dep, i, 1 );

       	for( j=0; (s1=localmem[0][j]) > -1; j++ )
		{
           	rootnode[s1] += (double)len[i][0] * eff[s1];
			eff[s1] *= 0.5;
/*
           	rootnode[s1] *= 0.5;
*/
			
		}
       	for( j=0; (s2=localmem[1][j]) > -1; j++ )
		{
           	rootnode[s2] +=  (double)len[i][1] * eff[s2];
			eff[s2] *= 0.5;
/*
           	rootnode[s2] *= 0.5;
*/
				
		}
	}
	for( i=0; i<nseq; i++ ) 
	{
#if 1 /* 97.9.29 */
		rootnode[i] += GETA3;
#endif
#if 0
		reporterr(       "### rootnode for %d = %f\n", i, rootnode[i] );
#endif
	}
#if 1
	total = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		total += rootnode[i];
	}
#else
	total = 1.0;
#endif
		
	for( i=0; i<nseq; i++ ) 
	{
		node[i] = rootnode[i] / total;
	}

#if 0
	reporterr(       "weight array in counteff_simple\n" );
	for( i=0; i<nseq; i++ )
		reporterr(       "%f\n", node[i] );
	printf( "\n" );
	exit( 1 );
#endif
	free( rootnode );
	free( eff );
	FreeIntMtx( localmem );
}

void counteff_simple_double_nostatic( int nseq, int ***topol, double **len, double *node )
{
    int i, j, s1, s2;
	double total;
	double *rootnode;
	double *eff;

	rootnode = AllocateDoubleVec( nseq );
	eff = AllocateDoubleVec( nseq );

	for( i=0; i<nseq; i++ ) // 2014/06/07, fu no eff wo sakeru.
	{
		if( len[i][0] < 0.0 ) 
		{
			reporterr( "WARNING: negative branch length %f, step %d-0\n", len[i][0], i );
			len[i][0] = 0.0;
		}
		if( len[i][1] < 0.0 ) 
		{
			reporterr( "WARNING: negative branch length %f, step %d-1\n", len[i][1], i );
			len[i][1] = 0.0;
		}
	}
#if DEBUG
	for( i=0; i<nseq-1; i++ )
	{
		reporterr( "\nstep %d, group 0\n", i );
		for( j=0; topol[i][0][j]!=-1; j++) reporterr( "%3d ", topol[i][0][j] );
		reporterr( "\n", i );
		reporterr( "step %d, group 1\n", i );
		for( j=0; topol[i][1][j]!=-1; j++) reporterr( "%3d ", topol[i][1][j] );
		reporterr( "\n", i );
		reporterr(       "len0 = %f\n", len[i][0] );
		reporterr(       "len1 = %f\n", len[i][1] );
	}
#endif
    for( i=0; i<nseq; i++ )
	{
		rootnode[i] = 0.0;
		eff[i] = 1.0;
/*
		rootnode[i] = 1.0;
*/
	}
   	for( i=0; i<nseq-1; i++ )
   	{
       	for( j=0; (s1=topol[i][0][j]) > -1; j++ )
		{
           	rootnode[s1] += (double)len[i][0] * eff[s1];
			eff[s1] *= 0.5;
/*
           	rootnode[s1] *= 0.5;
*/
			
		}
       	for( j=0; (s2=topol[i][1][j]) > -1; j++ )
		{
           	rootnode[s2] +=  (double)len[i][1] * eff[s2];
			eff[s2] *= 0.5;
/*
           	rootnode[s2] *= 0.5;
*/
				
		}
	}
	for( i=0; i<nseq; i++ ) 
	{
#if 1 /* 97.9.29 */
		rootnode[i] += GETA3;
#endif
#if 0
		reporterr(       "### rootnode for %d = %f\n", i, rootnode[i] );
#endif
	}
#if 1
	total = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		total += rootnode[i];
	}
#else
	total = 1.0;
#endif
		
	for( i=0; i<nseq; i++ ) 
	{
		node[i] = rootnode[i] / total;
	}

#if 0
	reporterr(       "weight array in counteff_simple\n" );
	for( i=0; i<nseq; i++ )
		reporterr(       "%f\n", node[i] );
	printf( "\n" );
	exit( 1 );
#endif
	free( rootnode );
	free( eff );
}

void counteff_simple( int nseq, int ***topol, double **len, double *node )
{
    int i, j, s1, s2;
	double total;
#if 0
	static double rootnode[M];
	static double eff[M];
#else
	double *rootnode;
	double *eff;
	rootnode = AllocateDoubleVec( nseq );
	eff = AllocateDoubleVec( nseq );
#endif

#if DEBUG
	for( i=0; i<nseq; i++ ){
		reporterr(       "len0 = %f\n", len[i][0] );
		reporterr(       "len1 = %f\n", len[i][1] );
	}
#endif
    for( i=0; i<nseq; i++ )
	{
		rootnode[i] = 0.0;
		eff[i] = 1.0;
/*
		rootnode[i] = 1.0;
*/
	}
   	for( i=0; i<nseq-1; i++ )
   	{
       	for( j=0; (s1=topol[i][0][j]) > -1; j++ )
		{
           	rootnode[s1] += len[i][0] * eff[s1];
			eff[s1] *= 0.5;
/*
           	rootnode[s1] *= 0.5;
*/
			
		}
       	for( j=0; (s2=topol[i][1][j]) > -1; j++ )
		{
           	rootnode[s2] +=  len[i][1] * eff[s2];
			eff[s2] *= 0.5;
/*
           	rootnode[s2] *= 0.5;
*/
				
		}
	}
	for( i=0; i<nseq; i++ ) 
	{
#if 1 /* 97.9.29 */
		rootnode[i] += GETA3;
#endif
#if 0
		reporterr(       "### rootnode for %d = %f\n", i, rootnode[i] );
#endif
	}
#if 1
	total = 0.0;
	for( i=0; i<nseq; i++ ) 
	{
		total += rootnode[i];
	}
#else
	total = 1.0;
#endif
		
	for( i=0; i<nseq; i++ ) 
	{
		node[i] = rootnode[i] / total;
	}

#if 0
	reporterr(       "weight array in counteff_simple\n" );
	for( i=0; i<nseq; i++ )
		reporterr(       "%f\n", node[i] );
	printf( "\n" );
	exit( 1 );
#endif
#if 1
	free( rootnode );
	free( eff );
#endif
}


void counteff( int nseq, int ***topol, double **len, double **node )
{
    int i, j, k, s1, s2;
	double rootnode[M];
	double eff[M];

	if( mix ) 
	{
		switch( weight )
		{
			case( 2 ): 
				weight = 3;
				break;
			case( 3 ): 
				weight = 2;
				break;
			default: 
				ErrorExit( "mix error" );
				break;
		}
	}

	if( weight == 2 )
	{
	    for( i=0; i<nseq; i++ ) rootnode[i] = 0;
    	for( i=0; i<nseq-2; i++ )
    	{
        	for( j=0; topol[i][0][j]>-1; j++ )
            	rootnode[topol[i][0][j]]++;
        	for( j=0; topol[i][1][j]>-1; j++ )
            	rootnode[topol[i][1][j]]++;
        	for( j=0; topol[i][0][j]>-1; j++ )
        	{
            	s1 = topol[i][0][j];
            	for( k=0; topol[i][1][k]>-1; k++ )
            	{
                	s2 = topol[i][1][k];
                	node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2] - 1;
            	}
        	}
    	}
    	for( j=0; topol[nseq-2][0][j]>-1; j++ )
    	{
        	s1 = topol[nseq-2][0][j];
        	for( k=0; topol[nseq-2][1][k]>-1; k++ )
        	{
            	s2 = topol[nseq-2][1][k];
            	node[MIN(s1,s2)][MAX(s1,s2)] = rootnode[s1] + rootnode[s2];
        	}
    	}
   		for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
   	   		node[i][j] = ipower( 0.5, (int)node[i][j] ) + geta2;
		for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ ) 
			node[j][i] = node[i][j];
	}

	if( weight == 3 )
	{
#if DEBUG
		for( i=0; i<nseq; i++ ){
			reporterr(       "len0 = %f\n", len[i][0] );
			reporterr(       "len1 = %f\n", len[i][1] );
		}
#endif
	    for( i=0; i<nseq; i++ )
		{
			rootnode[i] = 0.0;
			eff[i] = 1.0;
/*
			rootnode[i] = 1.0;
*/
		}
    	for( i=0; i<nseq-1; i++ )
    	{
        	for( j=0; (s1=topol[i][0][j]) > -1; j++ )
			{
   	        	rootnode[s1] += len[i][0] * eff[s1];
				eff[s1] *= 0.5;
/*
   	        	rootnode[s1] *= 0.5;
*/
				
			}
   	    	for( j=0; (s2=topol[i][1][j]) > -1; j++ )
			{
   	        	rootnode[s2] +=  len[i][1] * eff[s2];
				eff[s2] *= 0.5;
/*
   	        	rootnode[s2] *= 0.5;
*/
				
			}
		}
		for( i=0; i<nseq; i++ ) 
		{
#if 1 /* 97.9.29 */
			rootnode[i] += GETA3;
#endif
#if DEBUG
			reporterr(       "rootnode for %d = %f\n", i, rootnode[i] );
#endif
		}
		for( i=0; i<nseq; i++ ) 
		{
			for( j=0; j<nseq; j++ ) 
				if( j != i )
					node[i][j] = (double)rootnode[i] * rootnode[j];
				else node[i][i] = rootnode[i];
		}
	}

#if 0
	printf( "weight matrix in counteff\n" );
	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<nseq; j++ ) 
		{
			printf( "%f ", node[i][j] );
		}
		printf( "\n" );
	}
#endif
}

double score_calcp( char *seq1, char *seq2, int len )
{
	int k;
	unsigned char ms1, ms2;
	double tmpscore;
	int len2 = len - 2;

	tmpscore = 0.0;
	for( k=0; k<len; k++ )
	{
		ms1 = (unsigned char)seq1[k];
		ms2 = (unsigned char)seq2[k];
		if( ms1 == '-' && ms2 == '-' ) continue;
		tmpscore += (double)amino_dis[ms1][ms2];
	
		if( ms1 == (int)'-' ) 
		{
			tmpscore += (double)penalty;
			tmpscore += (double)amino_dis[ms1][ms2];
			while( (ms1=(unsigned char)seq1[++k]) == '-' )
				tmpscore += (double)amino_dis[ms1][ms2];
			k--;
			if( k >len2 ) break;
			continue;
		}
		if( ms2 == (int)'-' )
		{
			tmpscore += (double)penalty;
			tmpscore += (double)amino_dis[ms1][ms2];
			while( (ms2=(unsigned char)seq2[++k]) == '-' )
				tmpscore += (double)amino_dis[ms1][ms2];
			k--;
			if( k > len2 ) break;
			continue;
		}
	}
	return( tmpscore );
}

double score_calc1( char *seq1, char *seq2 )   /* method 1 */
{
	int k;
	double score = 0.0;
	int count = 0;
	int len = strlen( seq1 );

	for( k=0; k<len; k++ )
	{	
		if( seq1[k] != '-' && seq2[k] != '-' )
		{
			score += (double)amino_dis[(unsigned char)seq1[k]][(unsigned char)seq2[k]];
			count++;
		}
	}
	if( count ) score /= (double)count;
	else score = 1.0;
	return( score );
}

double substitution_nid( char *seq1, char *seq2 )
{
	int k;
	double s12;
	int len = strlen( seq1 );
	
	s12 = 0.0;
	for( k=0; k<len; k++ )
		if( seq1[k] != '-' && seq2[k] != '-' )
			s12 += ( seq1[k] == seq2[k] );

//	fprintf( stdout, "s12 = %f\n", s12 );
	return( s12 );
}

double substitution_score( char *seq1, char *seq2 )
{
	int k;
	double s12;
	int len = strlen( seq1 );
	
	s12 = 0.0;
	for( k=0; k<len; k++ )
		if( seq1[k] != '-' && seq2[k] != '-' )
			s12 += amino_dis[(unsigned char)seq1[k]][(unsigned char)seq2[k]];

//	fprintf( stdout, "s12 = %f\n", s12 );
	return( s12 );
}

double substitution_hosei( char *seq1, char *seq2 )   /* method 1 */
#if 0
{
	int k;
	double score = 0.0;
	int count = 0;
	int len = strlen( seq1 );

	for( k=0; k<len; k++ )
	{	
		if( seq1[k] != '-' && seq2[k] != '-' )
		{
			score += (double)( seq1[k] != seq2[k] );
			count++;
		}
	}
	if( count ) score /= (double)count;
	else score = 1.0;
	if( score < 0.95 ) score = - log( 1.0 - score );
	else score = 3.0;
	return( score );
}
#else
{
	int count = 0;
	double score;
	int iscore = 0;
	char s1, s2;

	while( (s1=*seq1++) )
	{
		s2 = *seq2++;
		if( s1 == '-' ) continue;
		if( s2 == '-' ) continue;
		iscore += ( s1 != s2 );
		count++;
	}
	if( count ) score = (double)iscore / count;
	else score = 1.0;
	if( score < 0.95 ) score = - log( 1.0 - score );
	else score = 3.0;
	return( score );
}
#endif

double substitution( char *seq1, char *seq2 )   /* method 1 */
{
	int k;
	double score = 0.0;
	int count = 0;
	int len = strlen( seq1 );

	for( k=0; k<len; k++ )
	{	
		if( seq1[k] != '-' && seq2[k] != '-' )
		{
			score += (double)( seq1[k] != seq2[k] );
			count++;
		}
	}
	if( count ) score /= (double)count;
	else score = 1.0;
	return( score );
}


void treeconstruction( char **seq, int nseq, int ***topol, double **len, double **eff )
{
    int i, j;

	if( weight > 1 )
	{
		if( utree == 0 )
		{
	    	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
   		 	{
/*
		       	 eff[i][j] = (double)score_calc1( seq[i], seq[j] );
*/
		       	 eff[i][j] = (double)substitution_hosei( seq[i], seq[j] );
 /*
				 reporterr(       "%f\n", eff[i][j] );
 */
   		 	}
/*
			reporterr(       "distance matrix\n" );
			for( i=0; i<nseq; i++ )
			{
				for( j=0; j<nseq; j++ ) 
				{
					reporterr(       "%f ", eff[i][j] );
				}
				reporterr(       "\n" );
			}
*/
/*
   			upg( nseq, eff, topol, len );
   			upg2( nseq, eff, topol, len );
*/
   			spg( nseq, eff, topol, len );
   			counteff( nseq, topol, len, eff );
		}
	}
	else
	{
		for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ ) 
			eff[i][j] = 1.0;
	}
/*
reporterr(       "weight matrix\n" );
for( i=0; i<nseq; i++ )
{
	for( j=0; j<nseq; j++ ) 
	{
		reporterr(       "%f ", eff[i][j] );
	}
	reporterr(       "\n" );
}
*/
}

double bscore_calc( char **seq, int s, double **eff )  /* algorithm B */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    long score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{
		double efficient = eff[i][j];

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 + !gb1  * !gc1 
                 * !gb2  *  gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
				 + gb1   * !gc1
				 * gb2   *  gc2      *BEFF

				 + gb1   *  gc1
				 * gb2   * !gc2      *BEFF
                 ;
			score += (long)cob * penalty * efficient;
			score += (long)amino_dis[(unsigned char)seq[i][k]][(unsigned char)seq[j][k]] * efficient;
			nglen += ( !gc1 * !gc2 );
		}
	}
	return( (double)score / nglen + 400.0 * !scoremtx );
}

void AllocateTmpSeqs( char ***mseq2pt, char **mseq1pt, int locnlenmax )
{
	*mseq2pt = AllocateCharMtx( njob, locnlenmax+1 );
	*mseq1pt = AllocateCharVec( locnlenmax+1 );
}

void FreeTmpSeqs( char **mseq2, char *mseq1 )
{
	FreeCharMtx( mseq2 );
	free( (char *)mseq1 );
}


void gappick0( char *aseq, char *seq )
{
	for( ; *seq != 0; seq++ )
	{
		if( *seq != '-' )
			*aseq++ = *seq;
	}
	*aseq = 0;

}

int isallgap( char *seq )
{
	for( ; *seq != 0; seq++ )
	{
		if( *seq != '-' )
			return( 0 );
	}
	return( 1 );
}

void gappick( int nseq, int s, char **aseq, char **mseq2, 
			  double **eff, double *effarr )
{
	int i, j, count, countjob, len, allgap;
	len = strlen( aseq[0] );
	for( i=0, count=0; i<len; i++ ) 
	{
		allgap = 1;
		for( j=0; j<nseq; j++ ) if( j != s ) allgap *= ( aseq[j][i] == '-' );
        if( allgap == 0 )
		{
			for( j=0, countjob=0; j<nseq; j++ ) 
			{
				if( j != s )
				{
					mseq2[countjob][count] = aseq[j][i];
					countjob++;
				}
			}
			count++;
		}
	}
	for( i=0; i<nseq-1; i++ ) mseq2[i][count] = 0;

	for( i=0, countjob=0; i<nseq; i++ ) 
	{
		if( i != s )
		{
			effarr[countjob] = eff[s][i];
			countjob++;
		}
	}
/*
fprintf( stdout, "effarr in gappick s = %d\n", s+1 );
for( i=0; i<countjob; i++ ) 
	fprintf( stdout, " %f", effarr[i] );
printf( "\n" );
*/
}

void commongappick_record( int nseq, char **seq, int *map )
{
	int i, j, count;
	int len = strlen( seq[0] );


	for( i=0, count=0; i<=len; i++ ) 
	{
	/*
		allgap = 1;
		for( j=0; j<nseq; j++ ) 
			allgap *= ( seq[j][i] == '-' );
		if( !allgap )
	*/
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j != nseq )
		{
			for( j=0; j<nseq; j++ )
			{
				seq[j][count] = seq[j][i];
			}
			map[count] = i;
			count++;
	 	}
	}
}


void commongappick( int nseq, char **seq )
{
	int i, j, count;
	int len = strlen( seq[0] );
#if 1

	int *mapfromnewtoold;

	mapfromnewtoold = calloc( len+1, sizeof( int ) );

	for( i=0, count=0; i<=len; i++ ) 
	{
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j != nseq )
		{
			mapfromnewtoold[count++] = i;
	 	}
	}
//	mapfromnewtoold[count] = -1; // iranai
	for( j=0; j<nseq; j++ )
	{
		for( i=0; i<count; i++ )
		{
			seq[j][i] = seq[j][mapfromnewtoold[i]];
		}
	}
	free( mapfromnewtoold );
#else
--------------------------

	int *mapfromoldtonew;
	int pos;

	mapfromoldtonew = calloc( len+1, sizeof( int ) );
	for( i=0; i<=len; i++ ) mapfromoldtonew[i] = -1;

	for( i=0, count=0; i<=len; i++ ) 
	{
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j != nseq )
		{
			mapfromoldtonew[i] = count;
			count++;
	 	}
	}
	for( j=0; j<nseq; j++ )
	{
		for( i=0; i<=len; i++ ) 
		{
			if( (pos=mapfromoldtonew[i]) != -1 )
				seq[j][pos] = seq[j][i];
		}
	}
	free( mapfromoldtonew );
--------------------------

	for( i=0, count=0; i<=len; i++ ) 
	{
	/*
		allgap = 1;
		for( j=0; j<nseq; j++ ) 
			allgap *= ( seq[j][i] == '-' );
		if( !allgap )
	*/
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j != nseq )
		{
			for( j=0; j<nseq; j++ )
			{
				seq[j][count] = seq[j][i];
			}
			count++;
	 	}
	}

#endif
}

#if 0
void commongaprecord( int nseq, char **seq, char *originallygapped )
{
	int i, j;
	int len = strlen( seq[0] );

	for( i=0; i<len; i++ ) 
	{
		for( j=0; j<nseq; j++ )
			if( seq[j][i] != '-' ) break;
		if( j == nseq )
			originallygapped[i] = '-';
		else
			originallygapped[i] = 'o';
	}
	originallygapped[len] = 0;
}
#endif
		
double score_calc0( char **seq, int s, double **eff, int ex )
{
	double tmp;

	if( scmtd == 4 ) tmp = score_calc4( seq, s, eff, ex );
	if( scmtd == 5 ) tmp = score_calc5( seq, s, eff, ex );
	else             tmp = score_calc5( seq, s, eff, ex );

	return( tmp );

}

/*
double score_m_1( char **seq, int ex, double **eff )
{
	int i, j, k;
	int len = strlen( seq[0] );
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
	double score;

	score = 0.0;
	nglen = 0;
	for( i=0; i<njob; i++ ) 
	{
		double efficient = eff[MIN(i,ex)][MAX(i,ex)];
		if( i == ex ) continue;

		gc1 = 0; 
		gc2 = 0;
		for( k=0; k<len; k++ ) 
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[ex][k] == '-' );
      
            cob = 
                   !gb1  *  gc1
                 * !gb2  * !gc2

                 + !gb1  * !gc1
                 * !gb2  *  gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
                 +  gb1  * !gc1
                 *  gb2  *  gc2      *BEFF

                 +  gb1  *  gc1
                 *  gb2  * !gc2      *BEFF
                 ;
			score += (double)cob * penalty * efficient;
			score += (double)amino_dis[seq[i][k]][seq[ex][k]] * efficient;
			*
			nglen += ( !gc1 * !gc2 );
			*
			if( !gc1 && !gc2 ) fprintf( stdout, "%f\n", score );
		}
	}
	return( (double)score / nglen + 400.0 * !scoremtx );
}
*/

#if 0
void sitescore( char **seq, double **eff, char sco1[], char sco2[], char sco3[] )
{
	int i, j, k;
	int len = strlen( seq[0] );
	double tmp;
	double count;
	int ch;
	double sco[N];

	for( i=0; i<len; i++ ) 
	{
		tmp = 0.0; count = 0;
		for( j=0; j<njob-1; j++ ) for( k=j+1; k<njob; k++ ) 
		{
		/*
			if( seq[j][i] != '-' && seq[k][i] != '-' )
		*/
			{
				tmp += amino_dis[seq[j][i]][seq[k][i]] + 400 * !scoremtx;
				count++; 
			}
		}
		if( count > 0.0 ) tmp /= count;
		else( tmp = 0.0 );
		ch = (int)( tmp/100.0 - 0.000001 );
		sprintf( sco1+i, "%c", ch+0x61 );
	}
	sco1[len] = 0;

    for( i=0; i<len; i++ ) 
    {
        tmp = 0.0; count = 0;
        for( j=0; j<njob-1; j++ ) for( k=j+1; k<njob; k++ ) 
        {
		/*
            if( seq[j][i] != '-' && seq[k][i] != '-' )
		*/
            {
                tmp += eff[j][k] * ( amino_dis[seq[j][i]][seq[k][i]] + 400 * !scoremtx );
                count += eff[j][k]; 
            }
        }
		if( count > 0.0 ) tmp /= count;
		else( tmp = 0.0 );
		tmp = ( tmp - 400 * !scoremtx ) * 2;
		if( tmp < 0 ) tmp = 0;
        ch = (int)( tmp/100.0 - 0.000001 );
        sprintf( sco2+i, "%c", ch+0x61 );
		sco[i] = tmp;
    }
    sco2[len] = 0;

	for( i=WIN; i<len-WIN; i++ )
	{
		tmp = 0.0;
		for( j=i-WIN; j<=i+WIN; j++ )
		{
			tmp += sco[j];
		}
		for( j=0; j<njob; j++ ) 
		{
			if( seq[j][i] == '-' )
			{
				tmp = 0.0;
				break;
			}
		}
		tmp /= WIN * 2 + 1;
		ch = (int)( tmp/100.0 - 0.0000001 );
		sprintf( sco3+i, "%c", ch+0x61 );
	}
	for( i=0; i<WIN; i++ ) sco3[i] = '-';
	for( i=len-WIN; i<len; i++ ) sco3[i] = '-';
	sco3[len] = 0;
}
#endif

void strins( char *str1, char *str2 )
{
	char *bk;
	int len1 = strlen( str1 );
	int len2 = strlen( str2 );

	bk = str2;
	str2 += len1+len2;
	str1 += len1-1;

	while( str2 >= bk+len1 ) { *str2 = *(str2-len1); str2--;} // by D.Mathog
	while( str2 >= bk ) { *str2-- = *str1--; }
}

int isaligned( int nseq, char **seq )
{
	int i;
	int len = strlen( seq[0] );
	for( i=1; i<nseq; i++ ) 
	{
		if( strlen( seq[i] ) != len ) return( 0 );
	}
	return( 1 );
}

double score_calc_for_score( int nseq, char **seq )
{
    int i, j, k, c;
    int len = strlen( seq[0] );
    double score;
    double tmpscore;
    char *mseq1, *mseq2;

    score = 0.0;
    for( i=0; i<nseq-1; i++ )
    {
        for( j=i+1; j<nseq; j++ )
        {
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0.0;
            c = 0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(unsigned char)mseq1[k]][(unsigned char)mseq2[k]];
                c++;
                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty - n_dis[0][24];
                    while( mseq1[++k] == '-' )
                        ;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty - n_dis[0][24];
                    while( mseq2[++k] == '-' )
                        ;
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            score += (double)tmpscore / (double)c;
#if DEBUG
			printf( "tmpscore in mltaln9.c = %f\n", tmpscore );
			printf( "tmpscore / c          = %f\n", tmpscore/(double)c );
#endif
        }
    }
	reporterr(       "raw score = %f\n", score );
	score /= (double)nseq * ( nseq-1.0 ) / 2.0;
	score += 400.0;
#if DEBUG
	printf( "score in mltaln9.c = %f\n", score );
#endif
    return( (double)score );
}

void doublencpy( double *vec1, double *vec2, int len )
{
	while( len-- )
		*vec1++ = *vec2++;
}

double score_calc_a( char **seq, int s, double **eff )  /* algorithm A+ */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    double score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{
		double efficient = eff[i][j];

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
				 +  gb1  * !gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 *  gb2  * !gc2
      
				 + !gb1  *  gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 * !gb2  *  gc2
                 ;
			score += 0.5 * (double)cob * penalty * efficient;
			score += (double)amino_dis[(unsigned char)seq[i][k]][(unsigned char)seq[j][k]] * (double)efficient;
			nglen += ( !gc1 * !gc2 );
		}
	}
	return( (double)score / nglen + 400.0 * !scoremtx );
}


double score_calc_s( char **seq, int s, double **eff )  /* algorithm S, not used */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    double score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{
		double efficient = eff[i][j];

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
#if 0
				 +  gb1  * !gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 *  gb2  * !gc2
      
				 + !gb1  *  gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 * !gb2  *  gc2
#endif
                 ;
			score += 0.5 * (double)cob * penalty * efficient;
			score += (double)amino_dis[(unsigned char)seq[i][k]][(int)seq[j][k]] * (double)efficient;
			nglen += ( !gc1 * !gc2 );
		}
	}
	return( (double)score / nglen + 400.0 );
}

double score_calc_for_score_s( int s, char **seq )  /* algorithm S */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    double score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
#if 0
				 +  gb1  * !gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 *  gb2  * !gc2
      
				 + !gb1  *  gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 * !gb2  *  gc2
#endif
                 ;
			score += 0.5 * (double)cob * penalty;
			score += (double)amino_dis[(int)seq[i][k]][(unsigned char)seq[j][k]];
			nglen += ( !gc1 * !gc2 );
		}
#if 0
		reporterr(       "i = %d, j=%d\n", i+1, j+1 );
		reporterr(       "score = %f\n", score );
#endif
	}
	return( (double)score / nglen + 400.0 );
}

double SSPscore___( int s, char **seq, int ex )  /* algorithm S */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    double score;

	score = 0;
	nglen = 0;
	i=ex; for( j=0; j<s; j++ )
	{

		if( j == ex ) continue;

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2 * 2.0 

                 +  gb1  * !gc1
                 * !gb2  *  gc2 * 2.0 
      
#if 0
				 +  gb1  * !gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 *  gb2  * !gc2
      
				 + !gb1  *  gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 * !gb2  *  gc2
#endif
                 ;
			score += 0.5 * (double)cob * penalty;
			score += (double)amino_dis[(unsigned char)seq[i][k]][(unsigned char)seq[j][k]];
			nglen += ( !gc1 * !gc2 ); /* tsukawanai */
		}
#if 0
		reporterr(       "i = %d, j=%d\n", i+1, j+1 );
		reporterr(       "score = %f\n", score );
#endif
	}
	return( (double)score );
}

double SSPscore( int s, char **seq )  /* algorithm S */
{
	int i, j, k;
	int gb1, gb2, gc1, gc2;
	int cob;
	int nglen;
    int len = strlen( seq[0] );
    double score;

	score = 0;
	nglen = 0;
	for( i=0; i<s-1; i++ ) for( j=i+1; j<s; j++ )
	{

		gc1 = 0;
		gc2 = 0;
		for( k=0; k<len; k++ )
		{
			gb1 = gc1;
			gb2 = gc2;

			gc1 = ( seq[i][k] == '-' );
			gc2 = ( seq[j][k] == '-' );
			
            cob = 
	               !gb1  *  gc1
 		         * !gb2  * !gc2

                 +  gb1  * !gc1 
                 * !gb2  * !gc2

	             + !gb1  * !gc1
 		         * !gb2  *  gc2

                 + !gb1  * !gc1 
                 *  gb2  * !gc2

                 + !gb1  *  gc1
                 *  gb2  * !gc2

                 +  gb1  * !gc1
                 * !gb2  *  gc2
      
#if 0
				 +  gb1  * !gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 *  gb2  * !gc2
      
				 + !gb1  *  gc1
				 *  gb2  *  gc2

				 +  gb1  *  gc1
				 * !gb2  *  gc2
#endif
                 ;
			score += 0.5 * (double)cob * penalty;
			score += (double)amino_dis[(unsigned char)seq[i][k]][(unsigned char)seq[j][k]];
			nglen += ( !gc1 * !gc2 ); /* tsukawanai */
		}
#if 0
		reporterr(       "i = %d, j=%d\n", i+1, j+1 );
		reporterr(       "score = %f\n", score );
#endif
	}
	return( (double)score );
}



double DSPscore( int s, char **seq )  /* method 3 deha nai */
{
    int i, j, k;
    double c;
    int len = strlen( seq[0] );
    double score;
    double tmpscore;
    char *mseq1, *mseq2;
#if DEBUG
	FILE *fp;
#endif

    score = 0.0;
    c = 0.0;

    for( i=0; i<s-1; i++ )
    {
        for( j=i+1; j<s; j++ )
        {
            mseq1 = seq[i];
            mseq2 = seq[j];
            tmpscore = 0.0;
            for( k=0; k<len; k++ )
            {
                if( mseq1[k] == '-' && mseq2[k] == '-' ) continue;
                tmpscore += amino_dis[(unsigned char)mseq1[k]][(unsigned char)mseq2[k]];

                if( mseq1[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq1[++k] == '-' )
                        tmpscore += amino_dis[(unsigned char)mseq1[k]][(unsigned char)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
                if( mseq2[k] == '-' )
                {
                    tmpscore += penalty;
                    while( mseq2[++k] == '-' )
                        tmpscore += amino_dis[(unsigned char)mseq1[k]][(unsigned char)mseq2[k]];
                    k--;
                    if( k > len-2 ) break;
                    continue;
                }
            }
            score += (double)tmpscore;
        }
    }

	return( score );
}


#define SEGMENTSIZE 150

int searchAnchors( int nseq, char **seq, Segment *seg )
{
	int i, j, k, kcyc;
	int status;
	double score;
	int value = 0;
	int len;
	int length;
	static double *stra = NULL;
	static int alloclen = 0;
	double cumscore;
	static double threshold;

	len = strlen( seq[0] );
	if( alloclen < len )
	{
		if( alloclen )
		{
			FreeDoubleVec( stra );
		}
		else
		{
			threshold = (int)divThreshold / 100.0 * 600.0 * divWinSize;
		}
		stra = AllocateDoubleVec( len );
		alloclen = len;
	}

	for( i=0; i<len; i++ )
	{
		stra[i] = 0.0;
		kcyc = nseq-1;
		for( k=0; k<kcyc; k++ ) for( j=k+1; j<nseq; j++ )
			stra[i] += n_dis[(int)amino_n[(unsigned char)seq[k][i]]][(int)amino_n[(unsigned char)seq[j][i]]];
		stra[i] /= (double)nseq * ( nseq-1 ) / 2;
	}

	(seg+0)->skipForeward = 0;
	(seg+1)->skipBackward = 0;
	status = 0;
	cumscore = 0.0;
	score = 0.0;
	length = 0; /* modified at 01/09/11 */
	for( j=0; j<divWinSize; j++ ) score += stra[j];
	for( i=1; i<len-divWinSize; i++ )
	{
		score = score - stra[i-1] + stra[i+divWinSize-1];
#if DEBUG
		reporterr(       "%d %f   ? %f", i, score, threshold );
		if( score > threshold ) reporterr(       "YES\n" );
		else                    reporterr(       "NO\n" );
#endif

		if( score > threshold )
		{
			if( !status )
			{
				status = 1;
				seg->start = i;
				length = 0;
				cumscore = 0.0;
			}
			length++;
			cumscore += score;
		}
		if( score <= threshold || length > SEGMENTSIZE )
		{
			if( status )
			{
				seg->end = i;
				seg->center = ( seg->start + seg->end + divWinSize ) / 2 ;
				seg->score = cumscore;
#if DEBUG
				reporterr(       "%d-%d length = %d\n", seg->start, seg->end, length );
#endif
				if( length > SEGMENTSIZE )
				{
					(seg+0)->skipForeward = 1;
					(seg+1)->skipBackward = 1;
				}
				else
				{
					(seg+0)->skipForeward = 0;
					(seg+1)->skipBackward = 0;
				}
				length = 0;
				cumscore = 0.0;
				status = 0;
				value++;
				seg++;
				if( value > MAXSEG - 3 ) ErrorExit( "TOO MANY SEGMENTS!");
			}
		}
	}
	if( status )
	{
		seg->end = i;
		seg->center = ( seg->start + seg->end + divWinSize ) / 2 ;
		seg->score = cumscore;
#if DEBUG
reporterr(       "%d-%d length = %d\n", seg->start, seg->end, length );
#endif
		value++;
	}
	return( value );
}

void dontcalcimportance_target( int nseq, double *eff, char **seq, LocalHom **localhom, int ntarget )
{
	int i, j;
	LocalHom *ptr;
	int *nogaplen;

	nogaplen = AllocateIntVec( nseq );

	for( i=0; i<nseq; i++ )
	{
		nogaplen[i] = seqlen( seq[i] );
//		reporterr(       "nogaplen[%d] = %d\n", i, nogaplen[i] );
	}

	for( i=0; i<ntarget; i++ )
	{
		for( j=0; j<nseq; j++ )
		{
			for( ptr=localhom[i]+j; ptr; ptr=ptr->next )
			{
//				reporterr(       "i,j=%d,%d,ptr=%p\n", i, j, ptr );
#if 1
				ptr->importance = ptr->opt / ptr->overlapaa;
//				ptr->fimportance = (double)ptr->importance;
#else
				ptr->importance = ptr->opt / MIN( nogaplen[i], nogaplen[j] );
#endif
			}
		}
	}
	free( nogaplen );
}
void dontcalcimportance( int nseq, double *eff, char **seq, LocalHom **localhom )
{
	int i, j;
	LocalHom *ptr;
	int *nogaplen;

	nogaplen = AllocateIntVec( nseq );

	for( i=0; i<nseq; i++ )
	{
		nogaplen[i] = seqlen( seq[i] );
//		reporterr(       "nogaplen[%d] = %d\n", i, nogaplen[i] );
	}

	for( i=0; i<nseq; i++ )
	{
		for( j=0; j<nseq; j++ )
		{
			for( ptr=localhom[i]+j; ptr; ptr=ptr->next )
			{
//				reporterr(       "i,j=%d,%d,ptr=%p\n", i, j, ptr );
#if 1
				ptr->importance = ptr->opt / ptr->overlapaa;
//				ptr->fimportance = (double)ptr->importance;
#else
				ptr->importance = ptr->opt / MIN( nogaplen[i], nogaplen[j] );
#endif
			}
		}
	}
	free( nogaplen );
}

void dontcalcimportance_firstone( int nseq, double *eff, char **seq, LocalHom **localhom )
{
	int i, j, nseq1;
	LocalHom *ptr;
#if 1
#else
	int *nogaplen;
	nogaplen = AllocateIntVec( nseq );
	for( i=0; i<nseq; i++ )
	{
		nogaplen[i] = seqlen( seq[i] );
//		reporterr(       "nogaplen[%d] = %d\n", i, nogaplen[i] );
	}
#endif

	nseq1 = nseq - 1;
	for( i=0; i<nseq1; i++ )
	{
		j=0;
		{
			for( ptr=localhom[i]+j; ptr; ptr=ptr->next )
			{
//				reporterr(       "i,j=%d,%d,ptr=%p\n", i, j, ptr );
#if 1
//				ptr->importance = ptr->opt / ptr->overlapaa;
				ptr->importance = ptr->opt * 0.5; // tekitou
//				ptr->fimportance = (double)ptr->importance;
//				reporterr(       "i=%d, j=%d, importance = %f, opt=%f\n", i, j, ptr->fimportance, ptr->opt  );
#else
				ptr->importance = ptr->opt / MIN( nogaplen[i], nogaplen[j] );
#endif
			}
		}
	}
#if 1
#else
	free( nogaplen );
#endif
}

void calcimportance_target( int nseq, int ntarget, double *eff, char **seq, LocalHom **localhom, int *targetmap, int *targetmapr )
{
	int i, j, pos, len, ti, tj;
	double *importance; // static -> local, 2012/02/25
	double tmpdouble;
	double *ieff, totaleff; // counteff_simple_double ni utsusu kamo
	int *nogaplen; // static -> local, 2012/02/25
	LocalHom *tmpptr;

	importance = AllocateDoubleVec( nlenmax );
	nogaplen = AllocateIntVec( nseq );
	ieff = AllocateDoubleVec( nseq );

	totaleff = 0.0;
	for( i=0; i<nseq; i++ )
	{
		nogaplen[i] = seqlen( seq[i] );
//		reporterr(       "nogaplen[] = %d\n", nogaplen[i] );
		if( nogaplen[i] == 0 ) ieff[i] = 0.0;
		else ieff[i] = eff[i];
		totaleff += ieff[i];
	}
	for( i=0; i<nseq; i++ ) ieff[i] /= totaleff;
	for( i=0; i<nseq; i++ ) printf(       "eff[%d] = %30.25f\n", i, ieff[i] );

#if 0
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		tmpptr = localhom[i]+j;
		reporterr(       "%d-%d\n", i, j );
		do
		{
			reporterr(       "reg1=%d-%d, reg2=%d-%d, opt=%f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt );
		} while( tmpptr=tmpptr->next );
	}
#endif


//	for( i=0; i<nseq; i++ )
	for( ti=0; ti<ntarget; ti++ )
	{
		i = targetmapr[ti];
//		reporterr(       "i = %d\n", i );
		for( pos=0; pos<nlenmax; pos++ )
			importance[pos] = 0.0;
		for( j=0; j<nseq; j++ )
		{
			if( i == j ) continue;
//			tmpptr = localhom[ti]+j;
			for( tmpptr = localhom[ti]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1 ) continue;
				for( pos=tmpptr->start1; pos<=tmpptr->end1; pos++ )
				{
#if 1
//					if( pos == 0 ) reporterr( "hit! i=%d, j=%d, pos=%d\n", i, j, pos );
					importance[pos] += ieff[j];
#else
					importance[pos] += ieff[j] * tmpptr->opt / MIN( nogaplen[i], nogaplen[j] );
					importance[pos] += ieff[j] * tmpptr->opt / tmpptr->overlapaa;
#endif
				}
			}
		}
#if 0
		reporterr(       "position specific importance of seq %d:\n", i );
		for( pos=0; pos<nlenmax; pos++ )
			reporterr(       "%d: %f\n", pos, importance[pos] );
		reporterr(       "\n" );
#endif
		for( j=0; j<nseq; j++ )
		{
//			reporterr(       "i=%d, j=%d\n", i, j );
			if( i == j ) continue;
			if( localhom[ti][j].opt == -1.0 ) continue;
#if 1
			for( tmpptr = localhom[ti]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1.0 ) continue;
				tmpdouble = 0.0;
				len = 0;
				for( pos=tmpptr->start1; pos<=tmpptr->end1; pos++ )
				{
					tmpdouble += importance[pos];
					len++;
				}

				tmpdouble /= (double)len;

				tmpptr->importance = tmpdouble * tmpptr->opt;
//				tmpptr->fimportance = (double)tmpptr->importance;
			}
#else
			tmpdouble = 0.0;
			len = 0;
			for( tmpptr = localhom[ti]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1.0 ) continue;
				for( pos=tmpptr->start1; pos<=tmpptr->end1; pos++ )
				{
					tmpdouble += importance[pos];
					len++;
				}
			}

			tmpdouble /= (double)len;

			for( tmpptr = localhom[ti]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1.0 ) continue;
				tmpptr->importance = tmpdouble * tmpptr->opt;
//				tmpptr->importance = tmpptr->opt / tmpptr->overlapaa; //$B$J$+$C$?$3$H$K$9$k(B
			}
#endif

//			reporterr(       "importance of match between %d - %d = %f\n", i, j, tmpdouble );
		}
	}

#if 0
	printf(       "before averaging:\n" );

	for( ti=0; ti<ntarget; ti++ ) for( j=0; j<nseq; j++ )
	{
		i = targetmapr[ti];
		if( i == j ) continue;
		printf(       "%d-%d\n", i, j );
		for( tmpptr = localhom[ti]+j; tmpptr; tmpptr=tmpptr->next )
		{
			printf(       "reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%30.25f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt / tmpptr->overlapaa, eff[i] * tmpptr->importance, tmpptr->opt );
		}
	}
#endif

#if 1
//	reporterr(       "average?\n" );
//	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
	for( ti=0; ti<ntarget; ti++ ) for( tj=ti+1; tj<ntarget; tj++ )
	{
		double imp;
		LocalHom *tmpptr1, *tmpptr2;

		i = targetmapr[ti];
		j = targetmapr[tj];
//		if( i == j ) continue;

//		reporterr(       "i=%d, j=%d\n", i, j );

		tmpptr1 = localhom[ti]+j; tmpptr2 = localhom[tj]+i;
		for( ; tmpptr1 && tmpptr2; tmpptr1 = tmpptr1->next, tmpptr2 = tmpptr2->next)
		{
			if( tmpptr1->opt == -1.0 || tmpptr2->opt == -1.0 ) 
			{
//				reporterr(       "WARNING: i=%d, j=%d, tmpptr1->opt=%f, tmpptr2->opt=%f\n", i, j, tmpptr1->opt, tmpptr2->opt );
				continue;
			}
//			reporterr(       "## importances = %f, %f\n", tmpptr1->importance, tmpptr2->importance );
			imp = 0.5 * ( tmpptr1->importance + tmpptr2->importance );
			tmpptr1->importance = tmpptr2->importance = imp;
//			tmpptr1->fimportance = tmpptr2->fimportance = (double)imp;

//			reporterr(       "## importance = %f\n", tmpptr1->importance );

		}

#if 0 // commented out, 2012/02/10
		if( ( tmpptr1 && !tmpptr2 ) || ( !tmpptr1 && tmpptr2 ) )
		{
			reporterr(       "ERROR: i=%d, j=%d\n", i, j );
			exit( 1 );
		}
#endif
	}

	for( ti=0; ti<ntarget; ti++ ) for( j=0; j<nseq; j++ )
	{
		double imp;
		LocalHom *tmpptr1;

		i = targetmapr[ti];
		if( i == j ) continue;
		if( targetmap[j] != -1 ) continue;

//		reporterr(       "i=%d, j=%d\n", i, j );

		tmpptr1 = localhom[ti]+j;
		for( ; tmpptr1; tmpptr1 = tmpptr1->next )
		{
			if( tmpptr1->opt == -1.0  ) 
			{
//				reporterr(       "WARNING: i=%d, j=%d, tmpptr1->opt=%f, tmpptr2->opt=%f\n", i, j, tmpptr1->opt, tmpptr2->opt );
				continue;
			}
//			reporterr(       "## importances = %f, %f\n", tmpptr1->importance, tmpptr2->importance );
			imp = 0.5 * ( tmpptr1->importance );
//			imp = 1.0 * ( tmpptr1->importance );
			tmpptr1->importance = imp;
//			tmpptr1->fimportance = (double)imp;

//			reporterr(       "## importance = %f\n", tmpptr1->importance );

		}

#if 0 // commented out, 2012/02/10
		if( ( tmpptr1 && !tmpptr2 ) || ( !tmpptr1 && tmpptr2 ) )
		{
			reporterr(       "ERROR: i=%d, j=%d\n", i, j );
			exit( 1 );
		}
#endif
	}
#endif
#if 0
	printf(       "after averaging:\n" );

	for( ti=0; ti<ntarget; ti++ ) for( j=0; j<nseq; j++ )
	{
		i = targetmapr[ti];
		if( i == j ) continue;
		printf(       "%d-%d\n", i, j );
		for( tmpptr = localhom[ti]+j; tmpptr; tmpptr=tmpptr->next )
		{
			if( tmpptr->end1 )
				printf(       "%d-%d, reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%f\n", i, j, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt / tmpptr->overlapaa, tmpptr->importance, tmpptr->opt );
		}
	}
//exit( 1 );
#endif
	free( importance );
	free( nogaplen );
	free( ieff );
}

void calcimportance_half( int nseq, double *eff, char **seq, LocalHom **localhom )
{
	int i, j, pos, len;
	double *importance; // static -> local, 2012/02/25
	double tmpdouble;
	double *ieff, totaleff; // counteff_simple_double ni utsusu kamo
	int *nogaplen; // static -> local, 2012/02/25
	LocalHom *tmpptr;

	importance = AllocateDoubleVec( nlenmax );
	nogaplen = AllocateIntVec( nseq );
	ieff = AllocateDoubleVec( nseq );

	totaleff = 0.0;
	for( i=0; i<nseq; i++ )
	{
		nogaplen[i] = seqlen( seq[i] );
//		reporterr(       "nogaplen[] = %d\n", nogaplen[i] );
		if( nogaplen[i] == 0 ) ieff[i] = 0.0;
		else ieff[i] = eff[i];
		totaleff += ieff[i];
	}
	for( i=0; i<nseq; i++ ) ieff[i] /= totaleff;
//	for( i=0; i<nseq; i++ ) reporterr(       "eff[%d] = %f\n", i, ieff[i] );

#if 0
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		tmpptr = localhom[i]+j;
		reporterr(       "%d-%d\n", i, j );
		do
		{
			reporterr(       "reg1=%d-%d, reg2=%d-%d, opt=%f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt );
		} while( tmpptr=tmpptr->next );
	}
#endif


	for( i=0; i<nseq; i++ )
	{
//		reporterr(       "i = %d\n", i );
		for( pos=0; pos<nlenmax; pos++ )
			importance[pos] = 0.0;
		for( j=0; j<nseq; j++ )
		{
			if( i == j ) continue;

			else if( i < j )
			{
				for( tmpptr = localhom[i]+j-i; tmpptr; tmpptr=tmpptr->next )
				{
					if( tmpptr->opt == -1 ) continue;
					for( pos=tmpptr->start1; pos<=tmpptr->end1; pos++ )
					{
#if 1
//						if( pos == 0 ) reporterr( "hit! i=%d, j=%d, pos=%d\n", i, j, pos );
						importance[pos] += ieff[j];
#else
						importance[pos] += ieff[j] * tmpptr->opt / MIN( nogaplen[i], nogaplen[j] );
						importance[pos] += ieff[j] * tmpptr->opt / tmpptr->overlapaa;
#endif
					}
				}
			}
			else
			{
				for( tmpptr = localhom[j]+i-j; tmpptr; tmpptr=tmpptr->next )
				{
					if( tmpptr->opt == -1 ) continue;
					for( pos=tmpptr->start2; pos<=tmpptr->end2; pos++ )
					{
#if 1
//						if( pos == 0 ) reporterr( "hit! i=%d, j=%d, pos=%d\n", i, j, pos );
						importance[pos] += ieff[j];
#else
						importance[pos] += ieff[j] * tmpptr->opt / MIN( nogaplen[i], nogaplen[j] );
						importance[pos] += ieff[j] * tmpptr->opt / tmpptr->overlapaa;
#endif
					}
				}
			}
		}
#if 0
		reporterr(       "position specific importance of seq %d:\n", i );
		for( pos=0; pos<nlenmax; pos++ )
			reporterr(       "%d: %f\n", pos, importance[pos] );
		reporterr(       "\n" );
#endif
		for( j=0; j<nseq; j++ )
		{
//			reporterr(       "i=%d, j=%d\n", i, j );
			if( i == j ) continue;

			else if( i < j )
			{
				if( localhom[i][j-i].opt == -1.0 ) continue;
	
				for( tmpptr = localhom[i]+j-i; tmpptr; tmpptr=tmpptr->next )
				{
					if( tmpptr->opt == -1.0 ) continue;
					tmpdouble = 0.0;
					len = 0;
					for( pos=tmpptr->start1; pos<=tmpptr->end1; pos++ )
					{
						tmpdouble += importance[pos];
						len++;
					}
	
					tmpdouble /= (double)len;
	
					tmpptr->importance = tmpdouble * tmpptr->opt;
//					tmpptr->fimportance = (double)tmpptr->importance;
				}
			}
			else
			{
				if( localhom[j][i-j].opt == -1.0 ) continue;
	
				for( tmpptr = localhom[j]+i-j; tmpptr; tmpptr=tmpptr->next )
				{
					if( tmpptr->opt == -1.0 ) continue;
					tmpdouble = 0.0;
					len = 0;
					for( pos=tmpptr->start2; pos<=tmpptr->end2; pos++ )
					{
						tmpdouble += importance[pos];
						len++;
					}
	
					tmpdouble /= (double)len;
	
					tmpptr->rimportance = tmpdouble * tmpptr->opt;
//					tmpptr->fimportance = (double)tmpptr->importance;
				}
			}

//			reporterr(       "importance of match between %d - %d = %f\n", i, j, tmpdouble );
		}
	}

#if 0
	printf(       "before averaging:\n" );

	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		if( i == j ) continue;

		else if( i < j ) 
		{
			printf(       "%d-%d\n", i, j );
			for( tmpptr = localhom[i]+j-i; tmpptr; tmpptr=tmpptr->next )
			{
				printf(       "reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt / tmpptr->overlapaa, eff[i] * tmpptr->importance, tmpptr->opt );
			}
		}
		else
		{
			printf(       "%d-%d\n", i, j );
			for( tmpptr = localhom[j]+i-j; tmpptr; tmpptr=tmpptr->next )
			{
				printf(       "reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%f\n", tmpptr->start2, tmpptr->end2, tmpptr->start1, tmpptr->end1, tmpptr->opt / tmpptr->overlapaa, eff[i] * tmpptr->rimportance, tmpptr->opt );
			}
		}
	}
#endif

#if 1
//	reporterr(       "average?\n" );
	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
	{
		double imp;
		LocalHom *tmpptr1;

//		reporterr(       "i=%d, j=%d\n", i, j );

		tmpptr1 = localhom[i]+j-i;
		for( ; tmpptr1; tmpptr1 = tmpptr1->next)
		{
			if( tmpptr1->opt == -1.0 )
			{
//				reporterr(       "WARNING: i=%d, j=%d, tmpptr1->opt=%f, tmpptr2->opt=%f\n", i, j, tmpptr1->opt, tmpptr2->opt );
				continue;
			}
//			reporterr(       "## importances = %f, %f\n", tmpptr1->importance, tmpptr2->importance );
			imp = 0.5 * ( tmpptr1->importance + tmpptr1->rimportance );
			tmpptr1->importance = tmpptr1->rimportance = imp;
//			tmpptr1->fimportance = tmpptr2->fimportance = (double)imp;

//			reporterr(       "## importance = %f\n", tmpptr1->importance );

		}

#if 0 // commented out, 2012/02/10
		if( ( tmpptr1 && !tmpptr2 ) || ( !tmpptr1 && tmpptr2 ) )
		{
			reporterr(       "ERROR: i=%d, j=%d\n", i, j );
			exit( 1 );
		}
#endif
	}
#endif
#if 0
	printf(       "after averaging:\n" );

	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		if( i < j ) for( tmpptr = localhom[i]+j-i; tmpptr; tmpptr=tmpptr->next )
		{
			if( tmpptr->end1 && tmpptr->start1 != -1 )
				printf(       "%d-%d, reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%f\n", i, j, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt / tmpptr->overlapaa, tmpptr->importance, tmpptr->opt );
		}
		else for( tmpptr = localhom[j]+i-j; tmpptr; tmpptr=tmpptr->next )
		{
			if( tmpptr->end2 && tmpptr->start2 != -1 )
				printf(       "%d-%d, reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%f\n", i, j, tmpptr->start2, tmpptr->end2, tmpptr->start1, tmpptr->end1, tmpptr->opt / tmpptr->overlapaa, tmpptr->importance, tmpptr->opt );
		}
	}
exit( 1 );
#endif
	free( importance );
	free( nogaplen );
	free( ieff );
}

void calcimportance( int nseq, double *eff, char **seq, LocalHom **localhom )
{
	int i, j, pos, len;
	double *importance; // static -> local, 2012/02/25
	double tmpdouble;
	double *ieff, totaleff; // counteff_simple_double ni utsusu kamo
	int *nogaplen; // static -> local, 2012/02/25
	LocalHom *tmpptr;

	importance = AllocateDoubleVec( nlenmax );
	nogaplen = AllocateIntVec( nseq );
	ieff = AllocateDoubleVec( nseq );

	totaleff = 0.0;
	for( i=0; i<nseq; i++ )
	{
		nogaplen[i] = seqlen( seq[i] );
//		reporterr(       "nogaplen[] = %d\n", nogaplen[i] );
		if( nogaplen[i] == 0 ) ieff[i] = 0.0;
		else ieff[i] = eff[i];
		totaleff += ieff[i];
	}
	for( i=0; i<nseq; i++ ) ieff[i] /= totaleff;
//	for( i=0; i<nseq; i++ ) reporterr(       "eff[%d] = %f\n", i, ieff[i] );

#if 0
	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		tmpptr = localhom[i]+j;
		reporterr(       "%d-%d\n", i, j );
		do
		{
			reporterr(       "reg1=%d-%d, reg2=%d-%d, opt=%f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt );
		} while( tmpptr=tmpptr->next );
	}
#endif


	for( i=0; i<nseq; i++ )
	{
//		reporterr(       "i = %d\n", i );
		for( pos=0; pos<nlenmax; pos++ )
			importance[pos] = 0.0;
		for( j=0; j<nseq; j++ )
		{
			if( i == j ) continue;
			tmpptr = localhom[i]+j;
			for( tmpptr = localhom[i]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1 ) continue;
				for( pos=tmpptr->start1; pos<=tmpptr->end1; pos++ )
				{
#if 1
//					if( pos == 0 ) reporterr( "hit! i=%d, j=%d, pos=%d\n", i, j, pos );
					importance[pos] += ieff[j];
#else
					importance[pos] += ieff[j] * tmpptr->opt / MIN( nogaplen[i], nogaplen[j] );
					importance[pos] += ieff[j] * tmpptr->opt / tmpptr->overlapaa;
#endif
				}
			}
		}
#if 0
		reporterr(       "position specific importance of seq %d:\n", i );
		for( pos=0; pos<nlenmax; pos++ )
			reporterr(       "%d: %f\n", pos, importance[pos] );
		reporterr(       "\n" );
#endif
		for( j=0; j<nseq; j++ )
		{
//			reporterr(       "i=%d, j=%d\n", i, j );
			if( i == j ) continue;
			if( localhom[i][j].opt == -1.0 ) continue;
#if 1
			for( tmpptr = localhom[i]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1.0 ) continue;
				tmpdouble = 0.0;
				len = 0;
				for( pos=tmpptr->start1; pos<=tmpptr->end1; pos++ )
				{
					tmpdouble += importance[pos];
					len++;
				}

				tmpdouble /= (double)len;

				tmpptr->importance = tmpdouble * tmpptr->opt;
//				tmpptr->fimportance = (double)tmpptr->importance;
			}
#else
			tmpdouble = 0.0;
			len = 0;
			for( tmpptr = localhom[i]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1.0 ) continue;
				for( pos=tmpptr->start1; pos<=tmpptr->end1; pos++ )
				{
					tmpdouble += importance[pos];
					len++;
				}
			}

			tmpdouble /= (double)len;

			for( tmpptr = localhom[i]+j; tmpptr; tmpptr=tmpptr->next )
			{
				if( tmpptr->opt == -1.0 ) continue;
				tmpptr->importance = tmpdouble * tmpptr->opt;
//				tmpptr->importance = tmpptr->opt / tmpptr->overlapaa; //$B$J$+$C$?$3$H$K$9$k(B
			}
#endif

//			reporterr(       "importance of match between %d - %d = %f\n", i, j, tmpdouble );
		}
	}

#if 0
	printf(       "before averaging:\n" );

	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		if( i == j ) continue;
		printf(       "%d-%d\n", i, j );
		for( tmpptr = localhom[i]+j; tmpptr; tmpptr=tmpptr->next )
		{
			printf(       "reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%f\n", tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt / tmpptr->overlapaa, eff[i] * tmpptr->importance, tmpptr->opt );
		}
	}
#endif

#if 1
//	reporterr(       "average?\n" );
	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
	{
		double imp;
		LocalHom *tmpptr1, *tmpptr2;

//		reporterr(       "i=%d, j=%d\n", i, j );

		tmpptr1 = localhom[i]+j; tmpptr2 = localhom[j]+i;
		for( ; tmpptr1 && tmpptr2; tmpptr1 = tmpptr1->next, tmpptr2 = tmpptr2->next)
		{
			if( tmpptr1->opt == -1.0 || tmpptr2->opt == -1.0 ) 
			{
//				reporterr(       "WARNING: i=%d, j=%d, tmpptr1->opt=%f, tmpptr2->opt=%f\n", i, j, tmpptr1->opt, tmpptr2->opt );
				continue;
			}
//			reporterr(       "## importances = %f, %f\n", tmpptr1->importance, tmpptr2->importance );
			imp = 0.5 * ( tmpptr1->importance + tmpptr2->importance );
			tmpptr1->importance = tmpptr2->importance = imp;
//			tmpptr1->fimportance = tmpptr2->fimportance = (double)imp;

//			reporterr(       "## importance = %f\n", tmpptr1->importance );

		}

#if 0 // commented out, 2012/02/10
		if( ( tmpptr1 && !tmpptr2 ) || ( !tmpptr1 && tmpptr2 ) )
		{
			reporterr(       "ERROR: i=%d, j=%d\n", i, j );
			exit( 1 );
		}
#endif
	}
#endif
#if 0
	printf(       "after averaging:\n" );

	for( i=0; i<nseq; i++ ) for( j=0; j<nseq; j++ )
	{
		for( tmpptr = localhom[i]+j; tmpptr; tmpptr=tmpptr->next )
		{
			if( tmpptr->end1 && tmpptr->start1 != -1 )
				printf(       "%d-%d, reg1=%d-%d, reg2=%d-%d, imp=%f -> %f opt=%f\n", i, j, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->opt / tmpptr->overlapaa, tmpptr->importance, tmpptr->opt );
		}
	}
exit( 1 );
#endif
	free( importance );
	free( nogaplen );
	free( ieff );
}




static void addlocalhom2_e( LocalHom *pt, LocalHom *lh, int sti, int stj, int eni, int enj, double opt, int overlp, int interm )
{
// dokka machigatteru
	if( pt != lh ) // susumeru
	{
		pt->next = (LocalHom *)calloc( 1, sizeof( LocalHom ) );
		pt = pt->next;
		pt->next = NULL;
		lh->last = pt;
	}
	else // sonomamatsukau
	{
		lh->last = pt;
	}
	lh->nokori++;
//	reporterr(       "in addlocalhom2_e, pt = %p, pt->next = %p, interm=%d, sti-eni-stj-enj=%d %d %d %d\n", pt, pt->next, interm, sti, eni, stj, enj );

	pt->start1 = sti;
	pt->start2 = stj;
	pt->end1 = eni;
	pt->end2 = enj;
	pt->opt = opt;
	pt->extended = interm;
	pt->overlapaa = overlp;
#if 0
	reporterr(       "i: %d-%d\n", sti, eni );
	reporterr(       "j: %d-%d\n", stj, enj );
	reporterr(       "opt=%f\n", opt );
	reporterr(       "overlp=%d\n", overlp );
#endif
}

void extendlocalhom2( int nseq, LocalHom **localhom, double **dist )
{
	int overlp, plim;
	int i, j, k;
	int pi, pj, pk, len;
	int status, sti, stj;
	int *ipt;
	int co;
	static int *ini = NULL;
	static int *inj = NULL;
	LocalHom *pt;

	sti = 0; // by D.Mathog, a guess
	stj = 0; // by D.Mathog, a guess

	if( ini == NULL )
	{
		ini = AllocateIntVec( nlenmax+1 );
		inj = AllocateIntVec( nlenmax+1 );
	}


	for( i=0; i<nseq-1; i++ )
	{
		for( j=i+1; j<nseq; j++ )
		{
#if 0
			for( k=0; k<nseq; k++ ) sai[k] = 0;
			numint = ncons;
			while( 1 )
			{
				k = (int)( rnd() * nseq );
				if( k == i || k == j ) continue; // mou yatta nomo habuita hoga ii 
				if( numint-- == 0 ) break;
				if( sai[k] ) continue;
				sai[k] = 1;
#else
			for( k=0; k<nseq; k++ )
			{
#endif
//				reporterr(       "i=%d, j=%d, k=%d, dists = %f,%f,%f thrinter=%f\n", i, j, k, dist[i][j], dist[MIN(i,k)][MAX(i,k)], dist[MIN(j,k)][MAX(j,k)], thrinter );
				if( k == i || k == j ) continue; // mou yatta nomo habuita hoga ii 
				if( dist[MIN(i,k)][MAX(i,k)] > dist[i][j] * thrinter || dist[MIN(j,k)][MAX(j,k)] > dist[i][j] * thrinter ) continue;
				ipt = ini; co = nlenmax+1;
				while( co-- ) *ipt++ = -1;
				ipt = inj; co = nlenmax+1;
				while( co-- ) *ipt++ = -1;
				overlp = 0;

				{
					for( pt=localhom[i]+k; pt; pt=pt->next )
		        	{
//						reporterr(       "i=%d,k=%d,st1:st2=%d:%d,pt=%p,extended=%p\n", i, k, pt->start1, pt->start2, pt, pt->extended );
						if( pt->opt == -1 )
						{
							reporterr(       "opt kainaide tbfast.c = %f\n", pt->opt );
						}
						if( pt->extended > -1 ) break;
						pi = pt->start1;
						pk = pt->start2;
						len = pt->end1 - pt->start1 + 1;
						ipt = ini + pk;
						while( len-- ) *ipt++ = pi++;
					}
				}

				{
					for( pt=localhom[j]+k; pt; pt=pt->next )
		        	{
						if( pt->opt == -1 )
						{
							reporterr(       "opt kainaide tbfast.c = %f\n", pt->opt );
						}
						if( pt->extended > -1 ) break;
						pj = pt->start1;
						pk = pt->start2;
						len = pt->end1 - pt->start1 + 1;
						ipt = inj + pk;
						while( len-- ) *ipt++ = pj++;
					}
				}
#if 0
				reporterr(       "i=%d,j=%d,k=%d\n", i, j, k );
				overlp = 0;
				for( pk = 0; pk < nlenmax; pk++ )
				{
					if( ini[pk] != -1 && inj[pk] != -1 ) overlp++;
					reporterr(       " %d", inj[pk] );
				}
				reporterr(       "\n" );

				reporterr(       "i=%d,j=%d,k=%d\n", i, j, k );
				overlp = 0;
				for( pk = 0; pk < nlenmax; pk++ )
				{
					if( ini[pk] != -1 && inj[pk] != -1 ) overlp++;
					reporterr(       " %d", ini[pk] );
				}
				reporterr(       "\n" );
#endif
				overlp = 0;
				plim = nlenmax+1;
				for( pk = 0; pk < plim; pk++ )
					if( ini[pk] != -1 && inj[pk] != -1 ) overlp++;


				status = 0;
				plim = nlenmax+1;
				for( pk=0; pk<plim; pk++ )
				{
//					reporterr(       "%d %d: %d-%d\n", i, j, ini[pk], inj[pk] );
					if( status )
					{
						if( ini[pk] == -1 || inj[pk] == -1 || ini[pk-1] != ini[pk] - 1 || inj[pk-1] != inj[pk] - 1 ) // saigonoshori
						{
							status = 0;
//							reporterr(       "end here!\n" );

							pt = localhom[i][j].last;
//							reporterr(       "in ex (ba), pt = %p, nokori=%d, i,j,k=%d,%d,%d\n", pt, localhom[i][j].nokori, i, j, k );
							addlocalhom2_e( pt, localhom[i]+j, sti, stj, ini[pk-1], inj[pk-1], MIN( localhom[i][k].opt, localhom[j][k].opt ) * 1.0, overlp, k );
//							reporterr(       "in ex, pt = %p, pt->next = %p, pt->next->next = %p\n", pt, pt->next, pt->next->next );

							pt = localhom[j][i].last;
//							reporterr(       "in ex (ba), pt = %p, pt->next = %p\n", pt, pt->next );
//							reporterr(       "in ex (ba), pt = %p, pt->next = %p, k=%d\n", pt, pt->next, k );
							addlocalhom2_e( pt, localhom[j]+i, stj, sti, inj[pk-1], ini[pk-1], MIN( localhom[i][k].opt, localhom[j][k].opt ) * 1.0, overlp, k );
//							reporterr(       "in ex, pt = %p, pt->next = %p, pt->next->next = %p\n", pt, pt->next, pt->next->next );
						}
					}
					if( !status ) // else deha arimasenn.
					{
						if( ini[pk] == -1 || inj[pk] == -1 ) continue;
						sti = ini[pk];
						stj = inj[pk];
//						reporterr(       "start here!\n" );
						status = 1;
					}
				}
//				if( status ) reporterr(       "end here\n" );

//				exit( 1 );
//					fprintf( hat3p, "%d %d %d %6.3f %d %d %d %d %p\n", i, j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->next ); 
			}
#if 0
			for( pt=localhomtable[i]+j; pt; pt=pt->next )
        	{
	            if( tmpptr->opt == -1.0 ) continue;
				fprintf( hat3p, "%d %d %d %6.3f %d %d %d %d %p\n", i, j, tmpptr->overlapaa, tmpptr->opt, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->next ); 
        	}
#endif
		}
	}
}

int makelocal( char *s1, char *s2, int thr )
{
	int start, maxstart, maxend;
	char *pt1, *pt2;
	double score;
	double maxscore;

	pt1 = s1;
	pt2 = s2;

	maxend = 0; // by D.Mathog, a guess

//	reporterr(       "thr = %d, \ns1 = %s\ns2 = %s\n", thr, s1, s2 );
	maxscore = 0.0;
	score = 0.0;
	start = 0;
	maxstart = 0;
	while( *pt1 )
	{
//		reporterr(       "*pt1 = %c*pt2 = %c\n", *pt1, *pt2 );
		if( *pt1 == '-' || *pt2 == '-' )
		{
//			reporterr(       "penalty = %d\n", penalty );
			score += penalty;
			while( *pt1 == '-' || *pt2 == '-' )
			{
				pt1++; pt2++;
			}
			continue;
		}

		score += ( amino_dis[(unsigned char)*pt1++][(unsigned char)*pt2++] - thr );
//		score += ( amino_dis[(int)*pt1++][(int)*pt2++] );
		if( score > maxscore ) 
		{
//			reporterr(       "score = %f\n", score );
			maxscore = score;
			maxstart = start;
//			reporterr(       "## max! maxstart = %d, start = %d\n", maxstart, start );
		}
		if( score < 0.0 )
		{
//			reporterr(       "## resetting, start = %d, maxstart = %d\n", start, maxstart );
			if( start == maxstart )
			{
				maxend = pt1 - s1;
//				reporterr(       "maxend = %d\n", maxend );
			}
			score = 0.0;
			start = pt1 - s1;
		}
	}
	if( start == maxstart )
		maxend = pt1 - s1 - 1;

//	reporterr(       "maxstart = %d, maxend = %d, maxscore = %f\n", maxstart, maxend, maxscore );
	s1[maxend+1] = 0;
	s2[maxend+1] = 0;
	return( maxstart );
}

void resetlocalhom( int nseq, LocalHom **lh )
{
	int i, j;
	LocalHom *pt;

	for( i=0; i<nseq-1; i++ ) for( j=i+1; j<nseq; j++ )
	{
		for( pt=lh[i]+j; pt; pt=pt->next )
			pt->opt = 1.0;
	}

}

void gapireru( char *res, char *ori, char *gt )
{
	char g;
	while( (g = *gt++) )
	{
		if( g == '-' )
		{
			*res++ = *newgapstr;
		}
		else
		{
			*res++ = *ori++;
		}
	}
	*res = 0;
}

void getkyokaigap( char *g, char **s, int pos, int n )
{
//	char *bk = g;
//	while( n-- ) *g++ = '-';
	while( n-- ) *g++ = (*s++)[pos];

//	reporterr(       "bk = %s\n", bk );
}

void new_OpeningGapCount( double *ogcp, int clus, char **seq, double *eff, int len, char *sgappat )
#if 0
{
	int i, j, gc, gb; 
	double feff;

	
	for( i=0; i<len+1; i++ ) ogcp[i] = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		gc = ( sgappat[j] == '-' );
		for( i=0; i<len; i++ ) 
		{
			gb = gc;
			gc = ( seq[j][i] == '-' );
			if( !gb *  gc ) ogcp[i] += feff;
		}
	}
}
#else
{
	int i, j, gc, gb; 
	double feff;
	double *fpt;
	char *spt;
	
	fpt = ogcp;
	i = len;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		spt = seq[j];
		fpt = ogcp;
		gc = ( sgappat[j] == '-' );
		i = len;
		while( i-- )
		{
			gb = gc;
			gc = ( *spt++ == '-' );
			{
				if( !gb *  gc ) *fpt += feff;
				fpt++;
			}
		}
	}
}
#endif
void new_OpeningGapCount_zure( double *ogcp, int clus, char **seq, double *eff, int len, char *sgappat, char *egappat )
#if 0
{
	int i, j, gc, gb; 
	double feff;

	
	for( i=0; i<len+1; i++ ) ogcp[i] = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		gc = ( sgappat[j] == '-' );
		for( i=0; i<len; i++ ) 
		{
			gb = gc;
			gc = ( seq[j][i] == '-' );
			if( !gb *  gc ) ogcp[i] += feff;
		}
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			if( !gb *  gc ) ogcp[i] += feff;
		}
	}
}
#else
{
	int i, j, gc, gb; 
	double feff;
	double *fpt;
	char *spt;
	
	fpt = ogcp;
	i = len+2;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		spt = seq[j];
		fpt = ogcp;
		gc = ( sgappat[j] == '-' );
		i = len;
		while( i-- )
		{
			gb = gc;
			gc = ( *spt++ == '-' );
			{
				if( !gb *  gc ) *fpt += feff;
				fpt++;
			}
		}
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			if( !gb *  gc ) *fpt += feff;
		}
	}
}
#endif

void new_FinalGapCount_zure( double *fgcp, int clus, char **seq, double *eff, int len, char *sgappat, char *egappat )
#if 0
{
	int i, j, gc, gb; 
	double feff;
	
	for( i=0; i<len+1; i++ ) fgcp[i] = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		gc = ( sgappat[j] == '-' );
		for( i=0; i<len; i++ ) 
		{
			gb = gc;
			gc = ( seq[j][i] == '-' );
			{
				if( gb * !gc ) fgcp[i] += feff;
			}
		}
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			{
				if( gb * !gc ) fgcp[len] += feff;
			}
		}
	}
}
#else
{
	int i, j, gc, gb; 
	double feff;
	double *fpt;
	char *spt;
	
	fpt = fgcp;
	i = len+2;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		fpt = fgcp;
		spt = seq[j];
		gc = ( sgappat[j] == '-' );
		i = len;
		while( i-- )
		{
			gb = gc;
			gc = ( *spt++ == '-' );
			{
				if( gb * !gc ) *fpt += feff;
				fpt++;
			}
		}
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			{
				if( gb * !gc ) *fpt += feff;
			}
		}
	}
}
#endif
void new_FinalGapCount( double *fgcp, int clus, char **seq, double *eff, int len, char *egappat )
#if 0
{
	int i, j, gc, gb; 
	double feff;
	
	for( i=0; i<len; i++ ) fgcp[i] = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		gc = ( seq[j][0] == '-' );
		for( i=1; i<len; i++ ) 
		{
			gb = gc;
			gc = ( seq[j][i] == '-' );
			{
				if( gb * !gc ) fgcp[i-1] += feff;
			}
		}
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			{
				if( gb * !gc ) fgcp[len-1] += feff;
			}
		}
	}
}
#else
{
	int i, j, gc, gb; 
	double feff;
	double *fpt;
	char *spt;
	
	fpt = fgcp;
	i = len;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		fpt = fgcp;
		spt = seq[j];
		gc = ( *spt == '-' );
		i = len;
		while( i-- )
		{
			gb = gc;
			gc = ( *++spt == '-' );
			{
				if( gb * !gc ) *fpt += feff;
				fpt++;
			}
		}
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			{
				if( gb * !gc ) *fpt += feff;
			}
		}
	}
}
#endif

void st_OpeningGapAdd( double *ogcp, int clus, char **seq, double *eff, int len )
{
	int i, j, gc, gb; 
	double *fpt;
	char *spt;
	int newmem = clus-1;
	double neweff = eff[newmem];
	double orieff = 1.0 - neweff;
	double feff;

//	fpt = ogcp;
//	i = len;
//	while( i-- ) *fpt++ = 0.0;
	
	j = clus-1;
//	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		spt = seq[j];
		fpt = ogcp;
		i = len;
		gc = 0;
		while( i-- )
		{
			gb = gc;
			gc = ( *spt++ == '-' );
			*fpt *= orieff;
			if( !gb *  gc ) *fpt += feff;
			fpt++;
		}
	}
	ogcp[len] = 0.0;

#if 0
	for( i=0; i<len; i++ )
		reporterr( "ogcp[%d]=%f\n", i, ogcp[i] );
	for( i=0; i<clus; i++ )
		reporterr( "%s\n", seq[i] );
	exit( 1 );
#endif
}

void st_OpeningGapCount( double *ogcp, int clus, char **seq, double *eff, int len )
{
	int i, j, gc, gb; 
	double feff;
	double *fpt;
	char *spt;
	
	fpt = ogcp;
	i = len;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		spt = seq[j];
		fpt = ogcp;
		gc = 0;
//		gc = 1;
		i = len;
		while( i-- )
		{
			gb = gc;
			gc = ( *spt++ == '-' );
			{
				if( !gb *  gc ) *fpt += feff;
				fpt++;
			}
		}
	}
	ogcp[len] = 0.0;
}

void st_FinalGapCount_zure( double *fgcp, int clus, char **seq, double *eff, int len )
{
	int i, j, gc, gb; 
	double feff;
	double *fpt;
	char *spt;
	
	fpt = fgcp;
	i = len+1;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		fpt = fgcp+1;
		spt = seq[j];
		gc = ( *spt == '-' );
		i = len;
//		for( i=1; i<len; i++ ) 
		while( i-- )
		{
			gb = gc;
			gc = ( *++spt == '-' );
			{
				if( gb * !gc ) *fpt += feff;
				fpt++;
			}
		}
		{
			gb = gc;
			gc = 0;
//			gc = 1;
			{
				if( gb * !gc ) *fpt += feff;
			}
		}
	}
}

void st_FinalGapAdd( double *fgcp, int clus, char **seq, double *eff, int len )
{
	int i, j, gc, gb; 
	double *fpt;
	char *spt;
	int newmem = clus-1;
	double neweff = eff[newmem];
	double orieff = 1.0 - neweff;
	double feff;
	
//	fpt = fgcp;
//	i = len;
//	while( i-- ) *fpt++ = 0.0;

	j = clus-1;
//	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		fpt = fgcp;
		spt = seq[j];
		gc = ( *spt == '-' );
		i = len;
//		for( i=1; i<len; i++ ) 
		while( i-- )
		{
			*fpt *= orieff;
			gb = gc;
			gc = ( *++spt == '-' );
			{
				if( gb * !gc ) *fpt += feff;
				fpt++;
			}
		}
		{
			*fpt *= orieff;
			gb = gc;
			gc = 0;
//			gc = 1;
			{
				if( gb * !gc ) *fpt += feff;
			}
		}
	}
}

void st_FinalGapCount( double *fgcp, int clus, char **seq, double *eff, int len )
{
	int i, j, gc, gb; 
	double feff;
	double *fpt;
	char *spt;
	
	fpt = fgcp;
	i = len;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		fpt = fgcp;
		spt = seq[j];
		gc = ( *spt == '-' );
		i = len;
//		for( i=1; i<len; i++ ) 
		while( i-- )
		{
			gb = gc;
			gc = ( *++spt == '-' );
			{
				if( gb * !gc ) *fpt += feff;
				fpt++;
			}
		}
		{
			gb = gc;
			gc = 0;
//			gc = 1;
			{
				if( gb * !gc ) *fpt += feff;
			}
		}
	}
}

void getGapPattern( double *fgcp, int clus, char **seq, double *eff, int len, char *xxx )
{
	int i, j, gc, gb; 
	double feff;
	double *fpt;
	char *spt;
	
	fpt = fgcp;
	i = len+1;
	while( i-- ) *fpt++ = 0.0;
	for( j=0; j<clus; j++ ) 
	{
		feff = (double)eff[j];
		fpt = fgcp;
		spt = seq[j];
		gc = ( *spt == '-' );
		i = len+1;
		while( i-- )
		{
			gb = gc;
			gc = ( *++spt == '-' );
			{
				if( gb * !gc ) *fpt += feff;
				fpt++;
			}
		}
#if 0
		{
			gb = gc;
			gc = ( egappat[j] == '-' );
			{
				if( gb * !gc ) *fpt += feff;
			}
		}
#endif
	}
	for( j=0; j<len; j++ )
	{
		reporterr(       "%d, %f\n", j, fgcp[j] );
	}
}

void getdigapfreq_st( double *freq, int clus, char **seq, double *eff, int len )
{
	int i, j;
	double feff;
	for( i=0; i<len+1; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		if( 0 && seq[i][0] == '-' ) // machigai kamo
			freq[0] += feff;
		for( j=1; j<len; j++ ) 
		{
			if( seq[i][j] == '-' && seq[i][j-1] == '-' )
				freq[j] += feff;
		}
		if( 0 && seq[i][len-1] == '-' )
			freq[len] += feff;
	}
//	reporterr(       "\ndigapf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getdiaminofreq_x( double *freq, int clus, char **seq, double *eff, int len )
{
	int i, j;
	double feff;
	for( i=0; i<len+2; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		if( seq[i][0] != '-' ) // tadashii
			freq[0] += feff;
		for( j=1; j<len; j++ ) 
		{
			if( seq[i][j] != '-' && seq[i][j-1] != '-' )
				freq[j] += feff;
		}
		if( 1 && seq[i][len-1] != '-' ) // xxx wo tsukawanaitoki [len-1] nomi
			freq[len] += feff;
	}
//	reporterr(       "\ndiaaf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getdiaminofreq_st( double *freq, int clus, char **seq, double *eff, int len )
{
	int i, j;
	double feff;
	for( i=0; i<len+1; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		if( seq[i][0] != '-' )
			freq[0] += feff;
		for( j=1; j<len; j++ ) 
		{
			if( seq[i][j] != '-' && seq[i][j-1] != '-' )
				freq[j] += feff;
		}
//		if( 1 && seq[i][len-1] != '-' ) // xxx wo tsukawanaitoki [len-1] nomi
			freq[len] += feff;
	}
//	reporterr(       "\ndiaaf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getdigapfreq_part( double *freq, int clus, char **seq, double *eff, int len, char *sgappat, char *egappat )
{
	int i, j;
	double feff;
	for( i=0; i<len+2; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
//		if( seq[i][0] == '-' )
		if( seq[i][0] == '-' && sgappat[i] == '-' )
			freq[0] += feff;
		for( j=1; j<len; j++ ) 
		{
			if( seq[i][j] == '-' && seq[i][j-1] == '-' )
				freq[j] += feff;
		}
//		if( seq[i][len] == '-' && seq[i][len-1] == '-' ) // xxx wo tsukawanaitoki arienai
		if( egappat[i] == '-' && seq[i][len-1] == '-' )
			freq[len] += feff;
	}
//	reporterr(       "\ndigapf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getdiaminofreq_part( double *freq, int clus, char **seq, double *eff, int len, char *sgappat, char *egappat )
{
	int i, j;
	double feff;
	for( i=0; i<len+2; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		if( seq[i][0] != '-' && sgappat[i] != '-' )
			freq[0] += feff;
		for( j=1; j<len; j++ ) 
		{
			if( seq[i][j] != '-' && seq[i][j-1] != '-' )
				freq[j] += feff;
		}
//		if( 1 && seq[i][len-1] != '-' ) // xxx wo tsukawanaitoki [len-1] nomi
		if( egappat[i] != '-'  && seq[i][len-1] != '-' ) // xxx wo tsukawanaitoki [len-1] nomi
			freq[len] += feff;
	}
//	reporterr(       "\ndiaaf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getgapfreq_zure_part( double *freq, int clus, char **seq, double *eff, int len, char *sgap )
{
	int i, j;
	double feff;
	for( i=0; i<len+2; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		if( sgap[i] == '-' )
			freq[0] += feff;
		for( j=0; j<len; j++ ) 
		{
			if( seq[i][j] == '-' )
				freq[j+1] += feff;
		}
//		if( egap[i] == '-' )
//			freq[len+1] += feff;
	}
//	reporterr(       "\ngapf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getgapfreq_zure( double *freq, int clus, char **seq, double *eff, int len )
{
	int i, j;
	double feff;
	for( i=0; i<len+1; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		for( j=0; j<len; j++ ) 
		{
			if( seq[i][j] == '-' )
				freq[j+1] += feff;
		}
	}
	freq[len+1] = 0.0;
//	reporterr(       "\ngapf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void getgapfreq( double *freq, int clus, char **seq, double *eff, int len )
{
	int i, j;
	double feff;
	for( i=0; i<len+1; i++ ) freq[i] = 0.0;
	for( i=0; i<clus; i++ )
	{
		feff = eff[i];
		for( j=0; j<len; j++ ) 
		{
			if( seq[i][j] == '-' )
				freq[j] += feff;
		}
	}
	freq[len] = 0.0;
//	reporterr(       "\ngapf = \n" );
//	for( i=0; i<len+1; i++ ) reporterr(       "%5.3f ", freq[i] );
}

void st_getGapPattern( Gappat **pat, int clus, char **seq, double *eff, int len )
{
	int i, j, k, gb, gc; 
	int known;
	double feff;
	Gappat **fpt;
	char *spt;
	int gaplen;

	fpt = pat;
	i = len+1;
	while( i-- ) 
	{
		if( *fpt ) free( *fpt );
		*fpt++ = NULL;
	}

	for( j=0; j<clus; j++ ) 
	{
//		reporterr(       "seq[%d] = %s\n", j, seq[j] );
		feff = (double)eff[j];

		fpt = pat;
		*fpt = NULL; // Falign.c kara yobareru tokiha chigau.
		spt = seq[j];
		gc = 0;
		gaplen = 0;

		for( i=0; i<len+1; i++ ) 
//		while( i-- )
		{
//			reporterr(       "i=%d, gaplen = %d\n", i, gaplen );
			gb = gc;
			gc = ( i != len && *spt++ == '-' );
			if( gc ) 
				gaplen++;
			else
			{
				if( gb && gaplen )
				{
					k = 1;
					known = 0;
					if( *fpt ) for( ; (*fpt)[k].len != -1; k++ )
					{
						if( (*fpt)[k].len == gaplen ) 
						{
//							reporterr(       "known\n" );
							known = 1;
							break;
						}
					}

					if( known == 0 )
					{
						*fpt = (Gappat *)realloc( *fpt, (k+3) *  sizeof( Gappat ) );  // mae1 (total), ato2 (len0), term
						if( !*fpt )
						{
							reporterr(       "Cannot allocate gappattern!'n" );
							reporterr(       "Use an approximate method, with the --mafft5 option.\n" );
							exit( 1 );
						}
						(*fpt)[k].freq = 0.0;
						(*fpt)[k].len = gaplen;
						(*fpt)[k+1].len = -1;
						(*fpt)[k+1].freq = 0.0; // iranai
//						reporterr(       "gaplen=%d, Unknown, %f\n", gaplen, (*fpt)[k].freq );
					}

//					reporterr(       "adding pos %d, len=%d, k=%d, freq=%f->", i, gaplen, k, (*fpt)[k].freq );
					(*fpt)[k].freq += feff;
//					reporterr(       "%f\n", (*fpt)[k].freq );
					gaplen = 0;
				}
			}
			fpt++;
		}
	}
#if 1
	for( j=0; j<len+1; j++ )
	{
		if( pat[j] )
		{
//			reporterr(       "j=%d\n", j );
//			for( i=1; pat[j][i].len!=-1; i++ )
//				reporterr(       "pos=%d, i=%d, len=%d, freq=%f\n", j, i, pat[j][i].len, pat[j][i].freq );

			pat[j][0].len = 0; // iminashi
			pat[j][0].freq = 0.0;
			for( i=1; pat[j][i].len!=-1;i++ )
			{
				pat[j][0].freq += pat[j][i].freq;
//				reporterr(       "totaling, i=%d, result = %f\n", i, pat[j][0].freq );
			}
//			reporterr(       "totaled, result = %f\n", pat[j][0].freq );

			pat[j][i].freq = 1.0 - pat[j][0].freq;
			pat[j][i].len = 0; // imiari
			pat[j][i+1].len = -1; 
		}
		else
		{
			pat[j] = (Gappat *)calloc( 3, sizeof( Gappat ) );
			pat[j][0].freq = 0.0;
			pat[j][0].len = 0; // iminashi

			pat[j][1].freq = 1.0 - pat[j][0].freq;
			pat[j][1].len = 0; // imiari
			pat[j][2].len = -1; 
		}
	}
#endif
}

static int minimum( int i1, int i2 )
{
	return MIN( i1, i2 );
}

static void commongappickpairfast( char *r1, char *r2, char *i1, char *i2, int *skip1, int *skip2 )
{
//	char *i1bk = i1;
	int skip, skipped1, skipped2;
//	int skip, skipped1, skipped2, scand1, scand2;
	skipped1 = skipped2 = 0;
//	reporterr("\n");
//	while( *i1 )
	while( 1 )
	{
//		fprintf( stderr, "i1 pos =%d\n", (int)(i1- i1bk) );
//		reporterr( "\nSkip cand %d-%d\n", *skip1-skipped1, *skip2-skipped2 );
#if 0
		scand1 = *skip1-skipped1;
		scand2 = *skip2-skipped2;
		skip = MIN( scand1, scand2 );
#else
		skip = minimum( *skip1-skipped1, *skip2-skipped2 );
#endif
//		reporterr( "Skip %d\n", skip );
		i1 += skip;
		i2 += skip;
		skipped1 += skip;
		skipped2 += skip;
//		fprintf( stderr, "i1 pos =%d, nlenmax=%d\n", (int)(i1- i1bk), nlenmax );
		if( !*i1 ) break;
//		reporterr( "%d, %c-%c\n", i1-i1bk, *i1, *i2 );
//		if( *i1 == '-' && *i2 == '-' )  // iranai?
//		{
//			reporterr( "Error in commongappickpairfast" );
//			exit( 1 );
//			i1++;
//			i2++;
//		}
		if( *i1 != '-' ) 
		{
			skipped1 = 0;
			skip1++;
		}
		else skipped1++;

		if( *i2 != '-' ) 
		{
			skipped2 = 0;
			skip2++;
		}
		else skipped2++;

		*r1++ = *i1++;
		*r2++ = *i2++;
	}
	*r1 = 0;
	*r2 = 0;
}

static void commongappickpair( char *r1, char *r2, char *i1, char *i2 )
{
//	strcpy( r1, i1 );
//	strcpy( r2, i2 );
//	return; // not SP
	while( *i1 )
	{
		if( *i1 == '-' && *i2 == '-' ) 
		{
			i1++;
			i2++;
		}
		else
		{
			*r1++ = *i1++;
			*r2++ = *i2++;
		}
	}
	*r1 = 0;
	*r2 = 0;
}

double naiveRpairscore( int n1, int n2, char **seq1, char **seq2, double *eff1, double *eff2, int penal )
{
//	return( 0 );
	int i, j;
	double val;
	double  valf;
	int  pv;
	double deff;
	char *p1, *p2, *p1p, *p2p;
	val = 0.0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		deff = eff1[i] * eff2[j];
//		reporterr(       "feff %d-%d = %f\n", i, j, feff );
//		reporterr(       "i1 = %s\n", seq1[i] );
//		reporterr(       "i2 = %s\n", seq2[j] );
//		reporterr(       "s1 = %s\n", s1 );
//		reporterr(       "s2 = %s\n", s2 );
//		reporterr(       "penal = %d\n", penal );

		valf = 0;
		p1 = seq1[i]; p2 = seq2[j];
		pv = 0;
		if( *p1 == '-' && *p2 != '-' )
			pv = penal;
		if( *p1 != '-' && *p2 == '-' )
			pv = penal;
//		if( pv ) reporterr(       "Penal!, %f, %d-%d, pos1,pos2=%d,%d\n", pv * deff * 0.5,  i, j, p1-seq1[i], p2-seq2[j] );
		p1p = p1; p2p = p2;
		valf += (double)amino_dis[(unsigned char)*p1++][(unsigned char)*p2++] + 0.5 * pv;
		while( *p1p )
		{
			pv = 0;
			if( *p1p != '-' && *p2p != '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					pv = penal;
				if( *p1 != '-' && *p2 == '-' )
					pv = penal;
				if( *p1 != '-' && *p2 != '-' )
					;
				if( *p1 == '-' && *p2 == '-' )
					;
			}
			if( *p1p == '-' && *p2p == '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					pv = penal;
//					;
				if( *p1 != '-' && *p2 == '-' )
					pv = penal;
//					;
				if( *p1 != '-' && *p2 != '-' )
					;
				if( *p1 == '-' && *p2 == '-' )
					;
			}
			if( *p1p != '-' && *p2p == '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					pv = penal * 2; // ??
//					;
				if( *p1 != '-' && *p2 == '-' )
					;
				if( *p1 != '-' && *p2 != '-' )
					pv = penal;
//					;
				if( *p1 == '-' && *p2 == '-' )
					pv = penal;
//					;
			}
			if( *p1p == '-' && *p2p != '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					;
				if( *p1 != '-' && *p2 == '-' )
					pv = penal * 2; // ??
//					;
				if( *p1 != '-' && *p2 != '-' )
					pv = penal;
//					;
				if( *p1 == '-' && *p2 == '-' )
					pv = penal;
//					;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
//			if( pv ) reporterr(       "Penal!, %f, %d-%d, pos1,pos2=%d,%d\n", pv * deff * 0.5,  i, j, p1-seq1[i], p2-seq2[j] );
			valf += amino_dis[(unsigned char)*p1++][(unsigned char)*p2++] + 0.5 * pv;
			p1p++; p2p++;
		}
//		reporterr(       "valf = %d\n", valf );
		val += deff * ( valf );
	}
	reporterr(       "val = %f\n", val );
	return( val );
//	exit( 1 );
}
double naiveQpairscore( int n1, int n2, char **seq1, char **seq2, double *eff1, double *eff2, int penal )
{
	int i, j;
	double val;
	double  valf;
	int  pv;
	double deff;
	char *p1, *p2, *p1p, *p2p;
	return( 0 );
	val = 0.0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		deff = eff1[i] * eff2[j];
//		reporterr(       "feff %d-%d = %f\n", i, j, feff );
//		reporterr(       "i1 = %s\n", seq1[i] );
//		reporterr(       "i2 = %s\n", seq2[j] );
//		reporterr(       "s1 = %s\n", s1 );
//		reporterr(       "s2 = %s\n", s2 );
//		reporterr(       "penal = %d\n", penal );

		valf = 0;
		p1 = seq1[i]; p2 = seq2[j];
		pv = 0;
		if( *p1 == '-' && *p2 != '-' )
			pv = penal;
		if( *p1 != '-' && *p2 == '-' )
			pv = penal;
//		if( pv ) reporterr(       "Penal!, %f, %d-%d, pos1,pos2=%d,%d\n", pv * deff * 0.5,  i, j, p1-seq1[i], p2-seq2[j] );
		p1p = p1; p2p = p2;
		valf += (double)amino_dis[(unsigned char)*p1++][(unsigned char)*p2++] + 0.5 * pv;
		while( *p1p )
		{
			pv = 0;
			if( *p1p != '-' && *p2p != '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					pv = penal;
				if( *p1 != '-' && *p2 == '-' )
					pv = penal;
				if( *p1 != '-' && *p2 != '-' )
					;
				if( *p1 == '-' && *p2 == '-' )
					;
			}
			if( *p1p == '-' && *p2p == '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
//					pv = penal;
					;
				if( *p1 != '-' && *p2 == '-' )
//					pv = penal;
					;
				if( *p1 != '-' && *p2 != '-' )
					;
				if( *p1 == '-' && *p2 == '-' )
					;
			}
			if( *p1p != '-' && *p2p == '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					pv = penal * 2; // ??
//					;
				if( *p1 != '-' && *p2 == '-' )
					;
				if( *p1 != '-' && *p2 != '-' )
					pv = penal;
//					;
				if( *p1 == '-' && *p2 == '-' )
//					pv = penal;
					;
			}
			if( *p1p == '-' && *p2p != '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					;
				if( *p1 != '-' && *p2 == '-' )
					pv = penal * 2; // ??
//					;
				if( *p1 != '-' && *p2 != '-' )
					pv = penal;
//					;
				if( *p1 == '-' && *p2 == '-' )
//					pv = penal;
					;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
//			if( pv ) reporterr(       "Penal!, %f, %d-%d, pos1,pos2=%d,%d\n", pv * deff * 0.5,  i, j, p1-seq1[i], p2-seq2[j] );
			valf += amino_dis[(unsigned char)*p1++][(unsigned char)*p2++] + 0.5 * pv;
			p1p++; p2p++;
		}
//		reporterr(       "valf = %d\n", valf );
		val += deff * ( valf );
	}
	reporterr(       "val = %f\n", val );
	return( val );
//	exit( 1 );
}
double naiveHpairscore( int n1, int n2, char **seq1, char **seq2, double *eff1, double *eff2, int penal )
{
	int i, j;
	double val;
	double  valf;
	int  pv;
//	double feff = 0.0; // by D.Mathog, a guess
	double deff;
	char *p1, *p2, *p1p, *p2p;
	val = 0.0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		deff = eff1[i] * eff2[j];
//		reporterr(       "i1 = %s\n", seq1[i] );
//		reporterr(       "i2 = %s\n", seq2[j] );
//		reporterr(       "s1 = %s\n", s1 );
//		reporterr(       "s2 = %s\n", s2 );
//		reporterr(       "penal = %d\n", penal );

		valf = 0;
		p1 = seq1[i]; p2 = seq2[j];
		pv = 0;
		if( *p1 == '-' && *p2 != '-' )
			pv = penal;
		if( *p1 != '-' && *p2 == '-' )
			pv = penal;
		if( pv ) reporterr(       "Penal!, %f, %d-%d, pos1,pos2=%d,%d\n", pv * deff * 0.5,  i, j, (int)(p1-seq1[i]), (int)(p2-seq2[j]) );
		p1p = p1; p2p = p2;
		valf += (double)amino_dis[(unsigned char)*p1++][(unsigned char)*p2++] + 0.5 * pv;
		while( *p1p )
		{
			pv = 0;
			if( *p1p != '-' && *p2p != '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					pv = penal;
				if( *p1 != '-' && *p2 == '-' )
					pv = penal;
				if( *p1 != '-' && *p2 != '-' )
					;
				if( *p1 == '-' && *p2 == '-' )
					;
			}
			if( *p1p == '-' && *p2p == '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
//					pv = penal;
					;
				if( *p1 != '-' && *p2 == '-' )
//					pv = penal;
					;
				if( *p1 != '-' && *p2 != '-' )
					;
				if( *p1 == '-' && *p2 == '-' )
					;
			}
			if( *p1p != '-' && *p2p == '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
//					pv = penal;
					;
				if( *p1 != '-' && *p2 == '-' )
					;
				if( *p1 != '-' && *p2 != '-' )
					pv = penal;
				if( *p1 == '-' && *p2 == '-' )
//					pv = penal;
					;
			}
			if( *p1p == '-' && *p2p != '-' )
			{
				if( *p1 == '-' && *p2 != '-' )
					;
				if( *p1 != '-' && *p2 == '-' )
//					pv = penal;
					;
				if( *p1 != '-' && *p2 != '-' )
					pv = penal;
				if( *p1 == '-' && *p2 == '-' )
//					pv = penal;
					;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
//			if( pv ) reporterr(       "Penal!, %f, %d-%d, pos1,pos2=%d,%d\n", pv * deff * 0.5,  i, j, p1-seq1[i], p2-seq2[j] );
			valf += amino_dis[(unsigned char)*p1++][(unsigned char)*p2++] + 0.5 * pv;
			p1p++; p2p++;
		}
//		reporterr(       "valf = %d\n", valf );
		val += deff * ( valf );
	}
	reporterr(       "val = %f\n", val );
	return( val );
//	exit( 1 );
}

double naivepairscorefast( char *seq1, char *seq2, int *skip1, int *skip2, int penal )
{
	double  vali;
	int len = strlen( seq1 );
	char *s1, *s2;
	char *p1, *p2;

	s1 = calloc( len+1, sizeof( char ) );
	s2 = calloc( len+1, sizeof( char ) );
	{
		vali = 0.0;
		commongappickpairfast( s1, s2, seq1, seq2, skip1, skip2 );
//		commongappickpair( s1, s2, seq1, seq2 );
//		printf(       "\n###s1 = %s\n", seq1 );
//		printf(       "###s2 = %s\n", seq2 );
//		printf(       "\n###i1 = %s\n", s1 );
//		printf(       "###i2 = %s\n", s2 );
//		printf( "allocated size, len+1 = %d\n", len+1 );
//		printf(       "###penal = %d\n", penal );

		p1 = s1; p2 = s2;
		while( *p1 )
		{
			if( *p1 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali += (double)penal;
//				while( *p1 == '-' || *p2 == '-' ) 
				while( *p1 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
			if( *p2 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali +=  (double)penal;
//				while( *p2 == '-' || *p1 == '-' ) 
				while( *p2 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
			vali += (double)amino_dis[(unsigned char)*p1++][(unsigned char)*p2++];
		}
	}
	free( s1 );
	free( s2 );
//	reporterr(       "###vali = %d\n", vali );
	return( vali );
}

double naivepairscore11_dynmtx( double **mtx, char *seq1, char *seq2, int penal )
{
	double  vali;
	int len = strlen( seq1 );
	char *s1, *s2, *p1, *p2;
	int c1, c2;


	s1 = calloc( len+1, sizeof( char ) );
	s2 = calloc( len+1, sizeof( char ) );
	{
		vali = 0.0;
		commongappickpair( s1, s2, seq1, seq2 );
//		reporterr(       "###i1 = %s\n", s1 );
//		reporterr(       "###i2 = %s\n", s2 );
//		reporterr(       "###penal = %d\n", penal );

		p1 = s1; p2 = s2;
		while( *p1 )
		{
			if( *p1 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali += (double)penal;
//				while( *p1 == '-' || *p2 == '-' ) 
				while( *p1 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
			if( *p2 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali +=  (double)penal;
//				while( *p2 == '-' || *p1 == '-' ) 
				while( *p2 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
			c1 = amino_n[(unsigned char)*p1++];
			c2 = amino_n[(unsigned char)*p2++];
			vali += (double)mtx[c1][c2];
		}
	}
	free( s1 );
	free( s2 );
//	reporterr(       "###vali = %d\n", vali );
	return( vali );
}

double naivepairscore11( char *seq1, char *seq2, int penal )
{
	double  vali;
	int len = strlen( seq1 );
	char *s1, *s2, *p1, *p2;

	s1 = calloc( len+1, sizeof( char ) );
	s2 = calloc( len+1, sizeof( char ) );
	{
		vali = 0.0;
		commongappickpair( s1, s2, seq1, seq2 );
//		reporterr(       "###i1 = %s\n", s1 );
//		reporterr(       "###i2 = %s\n", s2 );
//		reporterr(       "###penal = %d\n", penal );

		p1 = s1; p2 = s2;
		while( *p1 )
		{
			if( *p1 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali += (double)penal;
//				while( *p1 == '-' || *p2 == '-' ) 
				while( *p1 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
			if( *p2 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali +=  (double)penal;
//				while( *p2 == '-' || *p1 == '-' ) 
				while( *p2 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
			vali += (double)amino_dis[(unsigned char)*p1++][(unsigned char)*p2++];
		}
	}
	free( s1 );
	free( s2 );
//	reporterr(       "###vali = %d\n", vali );
	return( vali );
}

double naivepairscore( int n1, int n2, char **seq1, char **seq2, double *eff1, double *eff2, int penal )
{
//	return( 0.0 );
	int i, j;
	double val;
	int  vali;
	double feff;
	int len = strlen( seq1[0] );
	char *s1, *s2, *p1, *p2;
	s1 = calloc( len+1, sizeof( char ) );
	s2 = calloc( len+1, sizeof( char ) );
	val = 0.0;
	for( i=0; i<n1; i++ ) for( j=0; j<n2; j++ )
	{
		vali = 0;
		feff = eff1[i] * eff2[j];
//		reporterr(       "feff %d-%d = %f\n", i, j, feff );
		commongappickpair( s1, s2, seq1[i], seq2[j] );
//		reporterr(       "i1 = %s\n", seq1[i] );
//		reporterr(       "i2 = %s\n", seq2[j] );
//		reporterr(       "s1 = %s\n", s1 );
//		reporterr(       "s2 = %s\n", s2 );
//		reporterr(       "penal = %d\n", penal );

		p1 = s1; p2 = s2;
		while( *p1 )
		{
			if( *p1 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali += penal;
//				while( *p1 == '-' || *p2 == '-' ) 
				while( *p1 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
			if( *p2 == '-' )
			{
//				reporterr(       "Penal! %c-%c in %d-%d, %f\n", *(p1-1), *(p2-1), i, j, feff );
				vali +=  penal;
//				while( *p2 == '-' || *p1 == '-' ) 
				while( *p2 == '-' )  // SP
				{
					p1++;
					p2++;
				}
				continue;
			}
//			reporterr(       "adding %c-%c, %d\n", *p1, *p2, amino_dis[*p1][*p2] );
			vali += amino_dis[(unsigned char)*p1++][(unsigned char)*p2++];
		}
//		reporterr(       "vali = %d\n", vali );
		val += feff * vali;
	}
	free( s1 );
	free( s2 );
	reporterr(       "val = %f\n", val );
	return( val );
//	exit( 1 );
}

double plainscore( int nseq, char **s )
{
	int i, j, ilim;
	double v = 0.0;
	
	ilim = nseq-1;
	for( i=0; i<ilim; i++ ) for( j=i+1; j<nseq; j++ )
	{
		v += (double)naivepairscore11( s[i], s[j], penalty );
	}

	reporterr(       "penalty = %d\n", penalty );

	return( v );
}


int addonetip( int njobc, int ***topolc, double **lenc, double **iscorec, int ***topol, double **len, Treedep *dep, int treeout, Addtree *addtree, int iadd, char **name, int *alnleninnode, int *nogaplen, int noalign )
{
	int i, j, mem0, mem1, posinnew, m;
	int nstep;
	int norg;
	double minscore, minscoreo, eff0, eff1, addedlen, tmpmin;
	int nearest, nearesto;
	int repnorg;
	int *leaf2node;
	int *additionaltopol;
//	double (*clusterfuncpt[1])(double,double);
	Bchain *ac, *acpt, *acori, *acnext, *acprev;
	int neighbor;
	char *neighborlist;
	char *npt;
	int reflen, nearestnode, nogaplentoadd;
	int *topoldum0 = NULL;
	int *topoldum1 = NULL;
	int *topolo0;
	int *topolo1;
	int seqlengthcondition;
	double sueff1_double_local = 1.0 - sueff_global;
	double sueff05_double_local = sueff_global * 0.5;
//	char **tree; //static?
//	char *treetmp; //static?

//	for( i=0; i<njobc; i++ ) reporterr( "nogaplen of %d = %d\n", i+1, nogaplen[i] );
//exit( 1 );


//	treetmp = AllocateCharVec( njob*150 );
//	tree = AllocateCharMtx( njob, njob*150 );

//	sueff1_double = 1.0 - sueff_global;
//	sueff05_double = sueff_global * 0.5;
//	if ( treemethod == 'X' )
//		clusterfuncpt[0] = cluster_mix_double;
//	else if ( treemethod == 'E' )
//		clusterfuncpt[0] = cluster_average_double;
//	else if ( treemethod == 'q' )
//		clusterfuncpt[0] = cluster_minimum_double;
//	else
//	{
//		reporterr(       "Unknown treemethod, %c\n", treemethod );
//		exit( 1 );
//	}

	norg = njobc-1;
	nstep = njobc-2;

	additionaltopol = (int *)calloc( 2, sizeof( int ) );
	leaf2node= (int *)calloc( norg, sizeof( int ) );
	if( treeout )
	{
		neighborlist = calloc( norg * 30, sizeof( char ) );
	}
//	for( i=0; i<njobc; i++ ) sprintf( tree[i], "%d", i+1 );
	if( !leaf2node )
	{
		reporterr(       "Cannot allocate leaf2node.\n" );
		exit( 1 );
	}
	additionaltopol[0] = norg;
	additionaltopol[1] = -1;

	ac = (Bchain *)malloc( norg * sizeof( Bchain ) );
	for( i=0; i<norg; i++ )
	{
		ac[i].next = ac+i+1;
		ac[i].prev = ac+i-1;
		ac[i].pos = i;
	}
	ac[norg-1].next = NULL;


	acori = (Bchain *)malloc( 1 * sizeof( Bchain ) );
	acori->next = ac;
	acori->pos = -1;
	ac[0].prev = acori;


//	for( i=0; i<nstep; i++ )
//	{
//		reporterr(       "distfromtip = %f\n", dep[i].distfromtip );
//	}
//
//	for( i=0; i<norg; i++ )
//	{
//		reporterr(       "disttofrag(%d,%d) = %f\n", i, njobc-1, iscorec[i][norg-i] );
//	}


	minscore = 9999.9;
	nearest = -1;
	for( i=0; i<norg; i++ )
	{
		tmpmin = iscorec[i][norg-i];
		if( minscore > tmpmin )
		{
			minscore = tmpmin;
			nearest = i;
		}
	}
	nearesto = nearest;
	minscoreo = minscore;



//	for( i=0; i<njobc-1; i++ ) for( j=i+1; j<njobc; j++ )
//		reporterr(       "iscorec[%d][%d] = %f\n", i, j, iscorec[i][j-i] );
//	reporterr( "nearest = %d\n", nearest+1 );
//	reporterr( "nearesto = %d\n", nearesto+1 );

	posinnew = 0;
	repnorg = -1;
	nogaplentoadd = nogaplen[norg];



	for( i=0; i<norg; i++ ) leaf2node[i] = -1; 
	for( i=0; i<nstep; i++ )
	{
		mem0 = topol[i][0][0];
		mem1 = topol[i][1][0];
#if 0
		reporterr(       "\n\nstep %d (old) \n", i );

		reporterr( "group0 = \n" );
		for( j=0; topol[i][0][j]>-1; j++ ) 
		{
			reporterr( "%d ", topol[i][0][j]+1 );
		}
		reporterr( "\n" );
		reporterr( "len=%f\n", len[i][0] );
		reporterr( "group1 = \n" );
		for( j=0; topol[i][1][j]>-1; j++ ) 
		{
			reporterr( "%d ", topol[i][1][j]+1 );
		}
		reporterr( "\n" );
		reporterr( "len=%f\n", len[i][1] );
		
		reporterr(       "\n\n\nminscore = %f ? %f\n", minscore, dep[i].distfromtip*2 );
		reporterr(       "i = %d\n", i );
		if( leaf2node[nearest] == -1 )
		{
			reporterr(       "nogaplen[nearest] = %d\n", nogaplen[nearest] );
		}
		else
		{
			reporterr(       "alnleninnode[leaf2node[nearest]] = %d\n", alnleninnode[leaf2node[nearest]] );
			reporterr(       "leaf2node[nearest] = %d\n", leaf2node[nearest] );
		}
#endif
		nearestnode = leaf2node[nearest];
		if( nearestnode == -1 )
			reflen = nogaplen[nearest];
		else
			reflen = alnleninnode[nearestnode];
//			reflen = alnleninnode[i]; // BUG!!

		if( noalign ) seqlengthcondition =  1;
		else seqlengthcondition = ( nogaplentoadd <= reflen );

//seqlengthcondition = 1; // CHUUI
//seqlengthcondition = ( nogaplentoadd <= reflen ); // CHUUI

		if( repnorg == -1 && dep[i].distfromtip * 2 > minscore && seqlengthcondition )  // Keitouteki ichi ha fuseikaku.
//		if( repnorg == -1 && dep[i].distfromtip * 2 > minscore ) // Keitouteki ichi dake ga hitsuyouna baaiha kore wo tsukau.
		{
//			reporterr(       "INSERT HERE, %d-%d\n", nearest, norg );
//			reporterr(       "nearest = %d\n", nearest );
//			reporterr(       "\n\n\nminscore = %f\n", minscore );
//			reporterr(       "distfromtip *2 = %f\n", dep[i].distfromtip * 2 );
//			reporterr(       "nearest=%d, leaf2node[]=%d\n", nearest, leaf2node[nearest] );

			if( nearestnode == -1 )
			{
//				reporterr(       "INSERTING to 0!!!\n" );
//				reporterr(       "lastlength = %d\n", nogaplen[norg] );
//				reporterr(       "reflength = %d\n", nogaplen[nearest] );
				topolc[posinnew][0] = (int *)realloc( topolc[posinnew][0], ( 1 + 1 ) * sizeof( int ) );
				topolc[posinnew][0][0] = nearest;
				topolc[posinnew][0][1] = -1;

				addedlen = lenc[posinnew][0] = minscore / 2;

			}
			else
			{
//				reporterr(       "INSERTING to g, leaf2node = %d, cm=%d!!!\n", leaf2node[nearest], countmem(topol[leaf2node[nearest]][0] ) );
//				reporterr(       "alnleninnode[i] = %d\n", alnleninnode[i] );
//				reporterr(       "alnleninnode[leaf2node[nearest]] = %d\n", alnleninnode[leaf2node[nearest]] );

				topolc[posinnew][0] = (int *)realloc( topolc[posinnew][0], ( ( countmem( topol[nearestnode][0] ) + countmem( topol[nearestnode][1] ) + 1 ) * sizeof( int ) ) );
//				reporterr(       "leaf2node[%d] = %d\n", nearest, leaf2node[nearest] );
				intcpy( topolc[posinnew][0], topol[nearestnode][0] );
				intcat( topolc[posinnew][0], topol[nearestnode][1] );
//				addedlen = lenc[posinnew][0] = minscore / 2 - len[nearestnode][0]; // bug!!
				addedlen = lenc[posinnew][0] = dep[i].distfromtip - minscore / 2; // 2014/06/10
//				fprintf( stderr, "addedlen = %f, dep[i].distfromtip = %f, len[nearestnode][0] = %f, minscore/2 = %f, lenc[posinnew][0] = %f\n", addedlen, dep[i].distfromtip, len[nearestnode][0], minscore/2, lenc[posinnew][0] );

			}
			neighbor = lastmem( topolc[posinnew][0] );

			if( treeout )
			{
#if 0
				fp = fopen( "infile.tree", "a" ); // kyougou!!
				if( fp == 0 )
				{
					reporterr(       "File error!\n" );
					exit( 1 );
				}
				fprintf( fp, "\n" );
				fprintf( fp, "%8d: %s\n", norg+iadd+1, name[norg+iadd] );
				fprintf( fp, "          nearest sequence: %d\n", nearest + 1 );
				fprintf( fp, "          distance: %f\n", minscore );
				fprintf( fp, "          cousin: " );
				for( j=0; topolc[posinnew][0][j]!=-1; j++ )
					fprintf( fp, "%d ", topolc[posinnew][0][j]+1 );
				fprintf( fp, "\n" );
				fclose( fp );
#else
				addtree[iadd].nearest = nearesto;
				addtree[iadd].dist1 = minscoreo;
				addtree[iadd].dist2 = minscore;
				neighborlist[0] = 0;
				npt = neighborlist;
				for( j=0; topolc[posinnew][0][j]!=-1; j++ )
				{
					sprintf( npt, "%d ", topolc[posinnew][0][j]+1 );
					npt += strlen( npt );
				}
				addtree[iadd].neighbors = calloc( npt-neighborlist+1, sizeof( char ) );
				strcpy( addtree[iadd].neighbors, neighborlist );
#endif
			}

//			reporterr(       "INSERTING to 1!!!\n" );
			topolc[posinnew][1] = (int *)realloc( topolc[posinnew][1], ( 1 + 1 ) * sizeof( int ) );
			topolc[posinnew][1][0] = norg;
			topolc[posinnew][1][1] = -1;
			lenc[posinnew][1] = minscore / 2;

//			reporterr(       "STEP %d (newnew)\n", posinnew );
//			for( j=0; topolc[posinnew][0][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][0][j]+1 );
//			reporterr(       "\n len=%f\n", lenc[posinnew][0] );
//			for( j=0; topolc[posinnew][1][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][1][j]+1 );
//			reporterr(       "\n len=%f\n", lenc[posinnew][1] );

			repnorg = nearest;

//			reporterr(       "STEP %d\n", posinnew );
//			for( j=0; topolc[posinnew][0][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][0][j] );
//			reporterr(       "\n len=%f\n", lenc[i][0] );
//			for( j=0; topolc[posinnew][1][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][1][j] );
//			reporterr(       "\n len=%f\n", lenc[i][1] );

//			im = topolc[posinnew][0][0];
//			jm = topolc[posinnew][1][0];
//			sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], lenc[posinnew][0], tree[jm], lenc[posinnew][1] );
//			strcpy( tree[im], treetmp );

			posinnew++;
		}

//		reporterr(       "minscore = %f\n", minscore );
//		reporterr(       "distfromtip = %f\n", dep[i].distfromtip );
//		reporterr(       "Modify matrix, %d-%d\n", nearest, norg );
		eff0 = iscorec[mem0][norg-mem0];
		eff1 = iscorec[mem1][norg-mem1];

//		iscorec[mem0][norg-mem0] = (clusterfuncpt[0])( eff0, eff1 );
		iscorec[mem0][norg-mem0] =  MIN( eff0, eff1 ) * sueff1_double_local + ( eff0 + eff1 ) * sueff05_double_local; 
		iscorec[mem1][norg-mem1] = 9999.9; // sukoshi muda

		acprev = ac[mem1].prev; 
		acnext = ac[mem1].next; 
		acprev->next = acnext;
		if( acnext != NULL ) acnext->prev = acprev;

		if( ( nearest == mem1 || nearest == mem0 ) )
		{
			minscore = 9999.9;
//			for( j=0; j<norg; j++ ) // sukoshi muda
//			{
//				if( minscore > iscorec[j][norg-j] )
//				{
//					minscore = iscorec[j][norg-j];
//					nearest = j;
//				}
//			}
//			reporterr(       "searching on modified ac " );
			for( acpt=acori->next; acpt!=NULL; acpt=acpt->next ) // sukoshi muda
			{
//				reporterr(       "." );
				j = acpt->pos;
				tmpmin = iscorec[j][norg-j];
				if( minscore > tmpmin )
				{
					minscore = tmpmin;
					nearest = j;
				}
			}
//			reporterr(       "done\n" );
		}

//		reporterr(       "posinnew = %d\n", posinnew );


		if( topol[i][0][0] == repnorg )
		{
			topolc[posinnew][0] = (int *)realloc( topolc[posinnew][0], ( countmem( topol[i][0] ) + 2 ) * sizeof( int ) );
			intcpy( topolc[posinnew][0], topol[i][0] );
			intcat( topolc[posinnew][0], additionaltopol );
			lenc[posinnew][0] = len[i][0] - addedlen; // 2014/6/10
//			fprintf( stderr, "i=%d, dep[i].distfromtip=%f\n", i, dep[i].distfromtip );
//			fprintf( stderr, "addedlen=%f, len[i][0]=%f, lenc[][0]=%f\n", addedlen, len[i][0], lenc[posinnew][0] );
//			fprintf( stderr, "lenc[][1] = %f\n", lenc[posinnew][0] );
			addedlen = 0.0;
		}
		else
		{
			topolc[posinnew][0] = (int *)realloc( topolc[posinnew][0], ( countmem( topol[i][0] ) + 1 ) * sizeof( int ) );
			intcpy( topolc[posinnew][0], topol[i][0] );
			lenc[posinnew][0] = len[i][0];
		}

		if( topol[i][1][0] == repnorg )
		{
			topolc[posinnew][1] = (int *)realloc( topolc[posinnew][1], ( countmem( topol[i][1] ) + 2 ) * sizeof( int ) );
			intcpy( topolc[posinnew][1], topol[i][1] );
			intcat( topolc[posinnew][1], additionaltopol );
			lenc[posinnew][1] = len[i][1] - addedlen; // 2014/6/10
//			fprintf( stderr, "i=%d, dep[i].distfromtip=%f\n", i, dep[i].distfromtip );
//			fprintf( stderr, "addedlen=%f, len[i][1]=%f, lenc[][1]=%f\n", addedlen, len[i][1], lenc[posinnew][1] );
//			fprintf( stderr, "lenc[][1] = %f\n", lenc[posinnew][1] );
			addedlen = 0.0;

			repnorg = topolc[posinnew][0][0]; // juuyou
		}
		else
		{
			topolc[posinnew][1] = (int *)realloc( topolc[posinnew][1], ( countmem( topol[i][1] ) + 1 ) * sizeof( int ) );
			intcpy( topolc[posinnew][1], topol[i][1] );
			lenc[posinnew][1] = len[i][1];
		}

//		reporterr(       "\nSTEP %d (new)\n", posinnew );
//		for( j=0; topolc[posinnew][0][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][0][j]+1 );
//		reporterr(       "\n len=%f\n", lenc[posinnew][0] );
//		for( j=0; topolc[posinnew][1][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][1][j]+1 );
//		reporterr(       "\n len=%f\n", lenc[posinnew][1] );

//		reporterr("\ni=%d\n####### leaf2node[nearest]= %d\n", i, leaf2node[nearest] );

		for( j=0; (m=topol[i][0][j])!=-1; j++ ) leaf2node[m] = i;
		for( j=0; (m=topol[i][1][j])!=-1; j++ ) leaf2node[m] = i;

//		reporterr("####### leaf2node[nearest]= %d\n", leaf2node[nearest] );

//		im = topolc[posinnew][0][0];
//		jm = topolc[posinnew][1][0];
//		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], lenc[posinnew][0], tree[jm], lenc[posinnew][1] );
//		strcpy( tree[im], treetmp );
//
//		reporterr(       "%s\n", treetmp );

		posinnew++;
	}

	if( nstep )
	{
		i--;
		topolo0 = topol[i][0];
		topolo1 = topol[i][1];
	}
	else
	{
//		i = 0;
//		free( topol[i][0] );//?
//		free( topol[i][1] );//?
//		topol[i][0] = calloc( 2, sizeof( int ) );
//		topol[i][1] = calloc( 1, sizeof( int ) );
//		topol[i][0][0] = 0;
//		topol[i][0][1] = -1;
//		topol[i][1][0] = -1;

		topoldum0 = calloc( 2, sizeof( int ) );
		topoldum1 = calloc( 1, sizeof( int ) );
		topoldum0[0] = 0;
		topoldum0[1] = -1;
		topoldum1[0] = -1;

		topolo0 = topoldum0;
		topolo1 = topoldum1;
	}
	if( repnorg == -1 )
	{
//		topolc[posinnew][0] = (int *)realloc( topolc[posinnew][0], ( countmem( topol[i][0] ) + countmem( topol[i][1] ) + 1 ) * sizeof( int ) );
//		intcpy( topolc[posinnew][0], topol[i][0] );
//		intcat( topolc[posinnew][0], topol[i][1] );
		topolc[posinnew][0] = (int *)realloc( topolc[posinnew][0], ( countmem( topolo0 ) + countmem( topolo1 ) + 1 ) * sizeof( int ) );
		intcpy( topolc[posinnew][0], topolo0 );
		intcat( topolc[posinnew][0], topolo1 );
//		lenc[posinnew][0] = len[i][0] + len[i][1] - minscore / 2; // BUG!! 2014/06/07 ni hakken
		if( nstep )
			lenc[posinnew][0] = minscore / 2 - dep[nstep-1].distfromtip; // only when nstep>0, 2014/11/21
		else
			lenc[posinnew][0] = minscore / 2;

//		reporterr( "\ndep[nstep-1].distfromtip = %f\n", dep[nstep-1].distfromtip );
//		reporterr( "lenc[][0] = %f\n", lenc[posinnew][0] );

		topolc[posinnew][1] = (int *)realloc( topolc[posinnew][1],  2  * sizeof( int ) );
		intcpy( topolc[posinnew][1], additionaltopol );
		lenc[posinnew][1] = minscore / 2;

//		neighbor = lastmem( topolc[posinnew][0] );
		neighbor = norg-1;  // hakkirishita neighbor ga inai baai saigo ni hyouji

		if( treeout )
		{
#if 0
			fp = fopen( "infile.tree", "a" ); // kyougou!!
			if( fp == 0 )
			{
				reporterr(       "File error!\n" );
				exit( 1 );
			}
			fprintf( fp, "\n" );
			fprintf( fp, "%8d: %s\n", norg+iadd+1, name[norg+iadd] );
			fprintf( fp, "          nearest sequence: %d\n", nearest + 1 );
			fprintf( fp, "          cousin: " );
			for( j=0; topolc[posinnew][0][j]!=-1; j++ )
				fprintf( fp, "%d ", topolc[posinnew][0][j]+1 );
			fprintf( fp, "\n" );
			fclose( fp );
#else
			addtree[iadd].nearest = nearesto;
			addtree[iadd].dist1 = minscoreo;
			addtree[iadd].dist2 = minscore;
			neighborlist[0] = 0;
			npt = neighborlist;
			for( j=0; topolc[posinnew][0][j]!=-1; j++ )
			{
				sprintf( npt, "%d ", topolc[posinnew][0][j]+1 );
				npt += strlen( npt );
			}
			addtree[iadd].neighbors = calloc( npt-neighborlist+1, sizeof( char ) );
			strcpy( addtree[iadd].neighbors, neighborlist );
#endif
		}

//		reporterr(       "STEP %d\n", posinnew );
//		for( j=0; topolc[posinnew][0][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][0][j] );
//		reporterr(       "\n len=%f", lenc[posinnew][0] );
//		for( j=0; topolc[posinnew][1][j]!=-1; j++ ) reporterr(       " %d", topolc[posinnew][1][j] );
//		reporterr(       "\n len=%f\n", lenc[posinnew][1] );
	}

	if( topoldum0 ) free( topoldum0 );
	if( topoldum1 ) free( topoldum1 );
	free( leaf2node );
	free( additionaltopol );
	free( ac );
	free( acori );
	if( treeout ) free( neighborlist );

#if 0 //	create a newick tree for CHECK
	char **tree;
	char *treetmp;
	int im, jm;

	treetmp = AllocateCharVec( njob*150 );
	tree = AllocateCharMtx( njob, njob*150 );
	for( i=0; i<njobc; i++ ) sprintf( tree[i], "%d", i+1 );

	for( i=0; i<njobc-1; i++ )
	{
		reporterr( "\nSTEP %d\n", i );
		for( j=0; topolc[i][0][j]!=-1; j++ ) reporterr( " %d", topolc[i][0][j] );
		reporterr( "\n len=%f\n", lenc[i][0] );
		for( j=0; topolc[i][1][j]!=-1; j++ ) reporterr( " %d", topolc[i][1][j] );
		reporterr( "\n len=%f\n", lenc[i][1] );

		im = topolc[i][0][0];
		jm = topolc[i][1][0];
		sprintf( treetmp, "(%s:%7.5f,%s:%7.5f)", tree[im], lenc[i][0], tree[jm], lenc[i][1] );
		strcpy( tree[im], treetmp );

	}

	reporterr(       "%s\n", treetmp );
	FreeCharMtx( tree );
	free( treetmp );
#endif

	return( neighbor );
}

#if 0
int samemember( int *mem, int *cand )
{
	int i, j;

#if 0
	reporterr(       "mem = " );
	for( i=0; mem[i]>-1; i++ )	reporterr(       "%d ", mem[i] );
	reporterr(       "\n" );

	reporterr(       "cand = " );
	for( i=0; cand[i]>-1; i++ )	reporterr(       "%d ", cand[i] );
	reporterr(       "\n" );
#endif

	for( i=0, j=0; mem[i]>-1; )	
	{
		if( mem[i++] != cand[j++] ) return( 0 );
	}

	if( cand[j] == -1 )
	{
		return( 1 );
	}
	else
	{
		return( 0 );
	}
}
#else
int samemember( int *mem, int *cand )
{
	int i, j;
	int nm, nc;

	nm = 0; for( i=0; mem[i]>-1; i++ ) nm++;
	nc = 0; for( i=0; cand[i]>-1; i++ ) nc++;

	if( nm != nc ) return( 0 );

	for( i=0; mem[i]>-1; i++ )	
	{
		for( j=0; cand[j]>-1; j++ )
			if( mem[i] == cand[j] ) break;
		if( cand[j] == -1 ) return( 0 );
	}

	if( mem[i] == -1 )
	{
#if 0
		reporterr(       "mem = " );
		for( i=0; mem[i]>-1; i++ )	reporterr(       "%d ", mem[i] );
		reporterr(       "\n" );
	
		reporterr(       "cand = " );
		for( i=0; cand[i]>-1; i++ )	reporterr(       "%d ", cand[i] );
		reporterr(       "\n" );
#endif
		return( 1 );
	}
	else
	{
		return( 0 );
	}
}
#endif

int samemembern( int *mem, int *cand, int nc )
{
	int i, j;
	int nm;

	nm = 0; 
	for( i=0; mem[i]>-1; i++ ) 
	{
		nm++;
		if( nm > nc ) return( 0 );
	}

	if( nm != nc ) return( 0 );

	for( i=0; mem[i]>-1; i++ )	
	{
		for( j=0; j<nc; j++ )
			if( mem[i] == cand[j] ) break;
		if( j == nc ) return( 0 );
	}

	if( mem[i] == -1 )
	{
#if 0
		reporterr(       "mem = " );
		for( i=0; mem[i]>-1; i++ )	reporterr(       "%d ", mem[i] );
		reporterr(       "\n" );
	
		reporterr(       "cand = " );
		for( i=0; cand[i]>-1; i++ )	reporterr(       "%d ", cand[i] );
		reporterr(       "\n" );
#endif
		return( 1 );
	}
	else
	{
		return( 0 );
	}
}


int includemember( int *mem, int *cand ) // mem in cand 
{
	int i, j;

#if 0
	reporterr(       "mem = " );
	for( i=0; mem[i]>-1; i++ )	reporterr(       "%d ", mem[i] );
	reporterr(       "\n" );

	reporterr(       "cand = " );
	for( i=0; cand[i]>-1; i++ )	reporterr(       "%d ", cand[i] );
	reporterr(       "\n" );
#endif

	for( i=0; mem[i]>-1; i++ )
	{
		for( j=0; cand[j]>-1; j++ )
			if( mem[i] == cand[j] ) break;
		if( cand[j] == -1 ) return( 0 );
	}
//	reporterr(       "INCLUDED! mem[0]=%d\n", mem[0] );
	return( 1 );
}

int overlapmember( int *mem1, int *mem2 )
{
	int i, j;

	for( i=0; mem1[i]>-1; i++ )
		for( j=0; mem2[j]>-1; j++ )
			if( mem1[i] == mem2[j] ) return( 1 );
	return( 0 );
}
void gapcount( double *freq, char **seq, int nseq, double *eff, int lgth )
{
	int i, j;
	double fr;

//	for( i=0; i<lgth; i++ ) freq[i] = 0.0;
//	return;

	for( i=0; i<lgth; i++ )
	{
		fr = 0.0;
		for( j=0; j<nseq; j++ )
		{
			if( seq[j][i] == '-' ) fr += eff[j];
		}
		freq[i] = fr;
//		reporterr(       "freq[%d] = %f\n", i, freq[i] );
	}
//	reporterr(       "\n" );
	return;
}

void gapcountadd( double *freq, char **seq, int nseq, double *eff, int lgth )
{
	int i;
	int j = nseq-1;
	double newfr = eff[j];
	double orifr = 1.0 - newfr;

//	for( i=0; i<lgth; i++ ) freq[i] = 0.0;
//	return;
//	for( i=0; i<nseq; i++ )
//		reporterr( "%s\n", seq[i] );

	for( i=0; i<lgth; i++ )
	{
//		reporterr(       "freq[%d] = %f", i, freq[i] );
		freq[i] = 1.0 - freq[i]; // modosu
		freq[i] *= orifr;

		if( seq[j][i] == '-' ) freq[i] += newfr;
//		reporterr(       "->         %f\n", i, freq[i] );
	}
//	reporterr(       "\n" );
	return;
}
void gapcountf( double *freq, char **seq, int nseq, double *eff, int lgth )
{
	int i, j;
	double fr;

//	for( i=0; i<lgth; i++ ) freq[i] = 0.0;
//	return;

	for( i=0; i<lgth; i++ )
	{
		fr = 0.0;
		for( j=0; j<nseq; j++ )
		{
			if( seq[j][i] == '-' ) fr += eff[j];
		}
		freq[i] = fr;
//		reporterr(       "in gapcountf, freq[%d] = %f\n", i, freq[i] );
	}
//	reporterr(       "\n" );
	return;
}

void outgapcount( double *freq, int nseq, char *gappat, double *eff )
{
	int j;
	double fr;

	fr = 0.0;
	for( j=0; j<nseq; j++ )
	{
		if( gappat[j] == '-' ) fr += eff[j];
	}
	*freq = fr;
	return;
}

double dist2offset( double dist )
{
	double val = dist * 0.5 - specificityconsideration; // dist ha 0..2 dakara
//	double val = dist * 1.0 - specificityconsideration; // dist ha 0..2 dakara
	if( val > 0.0 ) val = 0.0;
	return val;
}

void makedynamicmtx( double **out, double **in, double offset )
{
	int i, j, ii, jj;
	double av;
 
	offset = dist2offset( offset * 2.0 ); // offset 0..1 -> 0..2

//	if( offset > 0.0 ) offset = 0.0;
//	reporterr(       "dynamic offset = %f\n", offset );

	for( i=0; i<nalphabets; i++ ) for( j=0; j<nalphabets; j++ )
	{
		out[i][j] = in[i][j];
	}
	if( offset == 0.0 ) return;

	for( i=0; i<nalphabets; i++ ) 
	{
		ii = (int)amino[i];
		if( ii == '-' ) continue; // text no toki arieru
		for( j=0; j<nalphabets; j++ )
		{
			jj = (int)amino[j];
			if( jj == '-' ) continue; // text no toki arieru
			out[i][j] = in[i][j] + offset * 600;
//			reporterr(       "%c-%c: %f\n", ii, jj, out[i][j] );
		}
	}

//	reporterr(       "offset = %f\n", offset );
//	reporterr(       "out[W][W] = %f\n", out[amino_n['W']][amino_n['W']] );
//	reporterr(       "out[A][A] = %f\n", out[amino_n['A']][amino_n['A']] );


	return;

// Taikaku youso no heikin ga 600 ni naruyouni re-scale.
// Hitaikaku youso ga ookiku narisugi.

	av = 0.0;
	for( i=0; i<nalphabets; i++ ) 
	{
		if( ii == '-' ) continue; // text no toki arieru
		av += out[i][i];
	}
	av /= (double)nalphabets;

	for( i=0; i<nalphabets; i++ ) 
	{
		if( amino[i] == '-' ) continue; // text no toki arieru
		for( j=0; j<nalphabets; j++ )
		{
			if( amino[j] == '-' ) continue; // text no toki arieru
			out[i][j] = out[i][j] * 600 / av;
			reporterr(       "%c-%c: %f\n", amino[i], amino[j], out[i][j] );
		}
	}
}
void FreeCommonIP()
{
	if( commonIP ) FreeIntMtx( commonIP ); 
	commonIP = NULL;
	commonAlloc1 = 0;
	commonAlloc2 = 0;
}

void makeskiptable( int n, int **skip, char **seq )
{
	char *nogapseq;
	int nogaplen, alnlen;
	int i, j, posinseq, gaplen;

	nogapseq = calloc( strlen( seq[0] )+1, sizeof( char ) );
	for( i=0; i<n; i++ )
	{
		gappick0( nogapseq, seq[i] );
		nogaplen = strlen( nogapseq );
		alnlen = strlen( seq[i] );
		skip[i] = calloc( nogaplen+1, sizeof( int ) );

//		reporterr( "%s\n", nogapseq );

		posinseq = 0;
		gaplen = 0;
		for( j=0; j<alnlen; j++ )
		{
			if( seq[i][j] == '-' )
			{
				skip[i][posinseq]++;
			}
			else
			{
				posinseq++;
			}
		}
//		for( j=0; j<nogaplen+1; j++ )
//			reporterr( "%d ", skip[i][j] );
//		reporterr( "\n" );
//		exit( 1 );
	}
	free( nogapseq );
}

int generatesubalignmentstable( int nseq, int ***tablept, int *nsubpt, int *maxmempt, int ***topol, double **len, double threshold )
{
	int i, j, rep0, rep1, nmem, mem;
	double distfromtip0, distfromtip1;
	double *distfromtip;
	reporterr( "\n\n\n" );

	*maxmempt = 0;
	*nsubpt = 0;

	distfromtip = calloc( nseq, sizeof( double ) );
	for( i=0; i<nseq-1; i++ )
	{
#if 0
		reporterr( "STEP %d\n", i );
		for( j=0; topol[i][0][j]!=-1; j++ )
			reporterr( "%3d ", topol[i][0][j] );
		reporterr( "\n" );
		reporterr( "len=%f\n", len[i][0] );
#endif

		rep0 = topol[i][0][0];
		distfromtip0 = distfromtip[rep0];
		distfromtip[rep0] += len[i][0];
//		reporterr( "distfromtip[%d] = %f->%f\n", rep0, distfromtip0, distfromtip[rep0] );


#if 0
		for( j=0; topol[i][1][j]!=-1; j++ )
			reporterr( "%3d ", topol[i][1][j] );
		reporterr( "\n" );
		reporterr( "len=%f\n", len[i][1] );
#endif

		rep1 = topol[i][1][0];
		distfromtip1 = distfromtip[rep1];
		distfromtip[rep1] += len[i][1];
//		reporterr( "distfromtip[%d] = %f->%f\n", rep1, distfromtip1, distfromtip[rep1] );

		if( topol[i][0][1] != -1 && distfromtip0 <= threshold && threshold < distfromtip[rep0] )
		{
//			reporterr( "HIT 0!\n" );
			*tablept = realloc( *tablept, sizeof( char * ) * (*nsubpt+2) );
			for( j=0, nmem=0; (mem=topol[i][0][j])!=-1; j++ ) 
				nmem++;
//			reporterr( "allocating %d\n", nmem+1 );
			(*tablept)[*nsubpt] = calloc( nmem+1, sizeof( int ) );
			(*tablept)[*nsubpt+1] = NULL;
			intcpy( (*tablept)[*nsubpt], topol[i][0] );
			if( *maxmempt < nmem ) *maxmempt = nmem;
			*nsubpt += 1;
		}

		if( topol[i][1][1] != -1 && distfromtip1 <= threshold && threshold < distfromtip[rep1] )
		{
//			reporterr( "HIT 1!\n" );
			*tablept = realloc( *tablept, sizeof( char * ) * (*nsubpt+2) );
			for( j=0, nmem=0; (mem=topol[i][1][j])!=-1; j++ )
				nmem++;
//			reporterr( "allocating %d\n", nmem+1 );
			(*tablept)[*nsubpt] = calloc( nmem+1, sizeof( int ) );
			(*tablept)[*nsubpt+1] = NULL;
			intcpy( (*tablept)[*nsubpt], topol[i][1] );
			if( *maxmempt < nmem ) *maxmempt = nmem;
			*nsubpt += 1;
		}

	}

	if( distfromtip[0] <= threshold ) 
	{
		free( distfromtip );
		return( 1 );
	}

	free( distfromtip );
	return( 0 );
}



double sumofpairsscore( int nseq, char **seq )
{
	double v = 0;
	int i, j;
	for( i=1; i<nseq; i++ )
	{
		for( j=0; j<i; j++ )
		{
			v += naivepairscore11( seq[i], seq[j], penalty ) / 600;
		}
	}
//	v /= ( (nseq-1) * nseq ) / 2;
	return( v );
}

int commonsextet_p( int *table, int *pointt )
{
	int value = 0;
	int tmp;
	int point;
	static TLS int *memo = NULL;
	static TLS int *ct = NULL;
	static TLS int *cp;

	if( table == NULL )
	{
		if( memo ) free( memo );
		if( ct ) free( ct );
		memo = NULL;
		ct = NULL;
		return( 0 );
	}

	if( *pointt == -1 )
		return( 0 );

	if( !memo )
	{
		memo = (int *)calloc( tsize, sizeof( int ) );
		if( !memo ) ErrorExit( "Cannot allocate memo\n" );
		ct = (int *)calloc( MIN( maxl, tsize )+1, sizeof( int ) ); // chuui!!
		if( !ct ) ErrorExit( "Cannot allocate ct\n" );
	}

	cp = ct;
	while( ( point = *pointt++ ) != END_OF_VEC )
	{
		tmp = memo[point]++;
		if( tmp < table[point] )
			value++;
		if( tmp == 0 ) *cp++ = point;
	}
	*cp = END_OF_VEC;
	
	cp =  ct;
	while( *cp != END_OF_VEC )
		memo[*cp++] = 0;

	return( value );
}

double distcompact_msa( char *seq1, char *seq2, int *skiptable1, int *skiptable2, int ss1, int ss2 ) // osoi!
{
	int bunbo = MIN( ss1, ss2 );
	double value;

//	reporterr( "msa-based dist\n" );
	if( bunbo == 0 )
		return( 2.0 );
	else
	{
		value = ( 1.0 - (double)naivepairscorefast( seq1, seq2, skiptable1, skiptable2, penalty_dist ) / bunbo ) * 2.0; // 2014/Aug/15 fast 
		if( value > 10 ) value = 10.0;  // 2015/Mar/17
		return( value );
	}
}

double distcompact( int len1, int len2, int *table1, int *point2, int ss1, int ss2 )
{
	double longer, shorter, lenfac, value;

	if( len1 > len2 )
	{
		longer=(double)len1;
		shorter=(double)len2;
	}
	else
	{
		longer=(double)len2;
		shorter=(double)len1;
	}
	lenfac = 1.0 / ( shorter / longer * lenfacd + lenfacb / ( longer + lenfacc ) + lenfaca );
//	reporterr( "lenfac=%f\n", lenfac );
//	reporterr( "commonsextet_p()=%d\n", commonsextet_p( table1, point2 ) );
//	reporterr( "ss1=%d, ss2=%d\n", ss1, ss2 );
//	reporterr( "val=%f\n", (1.0-(double)commonsextet_p( table1, point2 )/ss1) );

	if( ss1 == 0 || ss2 == 0 )
		return( 2.0 );

	value = ( 1.0 - (double)commonsextet_p( table1, point2 ) / MIN(ss1,ss2) ) * lenfac * 2.0;

	return( value ); // 2013/Oct/17 -> 2bai
}

static void movereg( char *seq1, char *seq2, LocalHom *tmpptr, int *start1pt, int *start2pt, int *end1pt, int *end2pt )
{
	char *pt;
	int tmpint;

	pt = seq1;
	tmpint = -1;
	while( *pt != 0 )
	{
		if( *pt++ != '-' ) tmpint++;
		if( tmpint == tmpptr->start1 ) break;
	}
	*start1pt = (int)( pt - seq1 ) - 1;
		
	if( tmpptr->start1 == tmpptr->end1 ) *end1pt = *start1pt;
	else
	{
		while( *pt != 0 )
		{
//			fprintf( stderr, "tmpint = %d, end1 = %d pos = %d\n", tmpint, tmpptr->end1, pt-seq1[i] );
			if( *pt++ != '-' ) tmpint++;
			if( tmpint == tmpptr->end1 ) break;
		}
		*end1pt = (int)( pt - seq1 ) - 1;
	}
		
	pt = seq2;
	tmpint = -1;
	while( *pt != 0 )
	{
		if( *pt++ != '-' ) tmpint++;
		if( tmpint == tmpptr->start2 ) break;
	}
	*start2pt = (int)( pt - seq2 ) - 1;
	if( tmpptr->start2 == tmpptr->end2 ) *end2pt = *start2pt;
	else
	{
		while( *pt != 0 )
		{
			if( *pt++ != '-' ) tmpint++;
			if( tmpint == tmpptr->end2 ) break;
		}
		*end2pt = (int)( pt - seq2 ) - 1;
	}
}

static void movereg_swap( char *seq1, char *seq2, LocalHom *tmpptr, int *start1pt, int *start2pt, int *end1pt, int *end2pt )
{
	char *pt;
	int tmpint;


	pt = seq1;
	tmpint = -1;
	while( *pt != 0 )
	{
		if( *pt++ != '-' ) tmpint++;
		if( tmpint == tmpptr->start2 ) break;
	}
	*start1pt = (int)( pt - seq1 ) - 1;
		
	if( tmpptr->start2 == tmpptr->end2 ) *end1pt = *start1pt;
	else
	{
		while( *pt != 0 )
		{
//			fprintf( stderr, "tmpint = %d, end1 = %d pos = %d\n", tmpint, tmpptr->end1, pt-seq1[i] );
			if( *pt++ != '-' ) tmpint++;
			if( tmpint == tmpptr->end2 ) break;
		}
		*end1pt = (int)( pt - seq1 ) - 1;
	}
		
	pt = seq2;
	tmpint = -1;
	while( *pt != 0 )
	{
		if( *pt++ != '-' ) tmpint++;
		if( tmpint == tmpptr->start1 ) break;
	}
	*start2pt = (int)( pt - seq2 ) - 1;
	if( tmpptr->start1 == tmpptr->end1 ) *end2pt = *start2pt;
	else
	{
		while( *pt != 0 )
		{
			if( *pt++ != '-' ) tmpint++;
			if( tmpint == tmpptr->end1 ) break;
		}
		*end2pt = (int)( pt - seq2 ) - 1;
	}
}

void fillimp( double **impmtx, double *imp, int clus1, int clus2, int lgth1, int lgth2, char **seq1, char **seq2, double *eff1, double *eff2, double *eff1_kozo, double *eff2_kozo, LocalHom ***localhom, char *swaplist, int forscore, int *orinum1, int *orinum2 )
{
	int i, j, k1, k2, start1, start2, end1, end2;
	double effij, effijx, effij_kozo; 
	char *pt1, *pt2;
	LocalHom *tmpptr;
	void (*movefunc)(char *, char *, LocalHom *, int *, int *, int *, int * );

#if 0
	fprintf( stderr, "eff1 in _init_strict = \n" );
	for( i=0; i<clus1; i++ )
		fprintf( stderr, "eff1[] = %f\n", eff1[i] );
	for( i=0; i<clus2; i++ )
		fprintf( stderr, "eff2[] = %f\n", eff2[i] );
#endif

	for( i=0; i<lgth1; i++ ) for( j=0; j<lgth2; j++ )
		impmtx[i][j] = 0.0;
	effijx = 1.0 * fastathreshold;
	for( i=0; i<clus1; i++ )
	{
		if( swaplist && swaplist[i] ) movefunc = movereg_swap;
		else movefunc = movereg;
		for( j=0; j<clus2; j++ )
		{

			if( swaplist == NULL && orinum1 && orinum2 ) // muda. 
			{
				if( orinum1[i]>orinum2[j] )
					movefunc = movereg_swap;
				else
					movefunc = movereg;
			}

//			effij = eff1[i] * eff2[j] * effijx;
			effij = eff1[i] * eff2[j] * effijx;
			effij_kozo = eff1_kozo[i] * eff2_kozo[j] * effijx;
			tmpptr = localhom[i][j];
			while( tmpptr )
			{
//				fprintf( stderr, "start1 = %d\n", tmpptr->start1 );
//				fprintf( stderr, "end1   = %d\n", tmpptr->end1   );
//				fprintf( stderr, "i = %d, seq1 = \n%s\n", i, seq1[i] );
//				fprintf( stderr, "j = %d, seq2 = \n%s\n", j, seq2[j] );

				movefunc( seq1[i], seq2[j], tmpptr, &start1, &start2, &end1, &end2 );


//				fprintf( stderr, "start1 = %d (%c), end1 = %d (%c), start2 = %d (%c), end2 = %d (%c)\n", start1, seq1[i][start1], end1, seq1[i][end1], start2, seq2[j][start2], end2, seq2[j][end2] );
//				fprintf( stderr, "step 0\n" );
				if( end1 - start1 != end2 - start2 )
				{
//					fprintf( stderr, "CHUUI!!, start1 = %d, end1 = %d, start2 = %d, end2 = %d\n", start1, end1, start2, end2 );
				}

				k1 = start1; k2 = start2;
				pt1 = seq1[i] + k1;
				pt2 = seq2[j] + k2;
				while( *pt1 && *pt2 )
				{
					if( *pt1 != '-' && *pt2 != '-' )
					{
// ½Å¤ß¤òÆó½Å¤Ë¤«¤±¤Ê¤¤¤è¤¦¤ËÃí°Õ¤·¤Æ²¼¤µ¤¤¡£
//						impmtx[k1][k2] += tmpptr->wimportance * fastathreshold;
//						impmtx[k1][k2] += tmpptr->importance * effij;
//						impmtx[k1][k2] += tmpptr->fimportance * effij;
						if( tmpptr->korh == 'k' )
							impmtx[k1][k2] += tmpptr->importance * effij_kozo;
						else
							impmtx[k1][k2] += tmpptr->importance * effij;
//						fprintf( stderr, "k1=%d, k2=%d, impalloclen=%d\n", k1, k2, impalloclen );
//						fprintf( stderr, "mark, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
						k1++; k2++;
						pt1++; pt2++;
					}
					else if( *pt1 != '-' && *pt2 == '-' )
					{
//						fprintf( stderr, "skip, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
						k2++; pt2++;
					}
					else if( *pt1 == '-' && *pt2 != '-' )
					{
//						fprintf( stderr, "skip, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
						k1++; pt1++;
					}
					else if( *pt1 == '-' && *pt2 == '-' )
					{
//						fprintf( stderr, "skip, %d (%c) - %d (%c) \n", k1, *pt1, k2, *pt2 );
						k1++; pt1++;
						k2++; pt2++;
					}
					if( k1 > end1 || k2 > end2 ) break;
				}
				tmpptr = tmpptr->next;
			}
		}
	}
#if 0
	printf( "orinum1=%d, orinum2=%d\n", *orinum1, *orinum2 );
	if( *orinum1 == 0 )
	{
		fprintf( stdout, "impmtx = \n" );
		for( k2=0; k2<lgth2; k2++ )
			fprintf( stdout, "%6.3f ", (double)k2 );
		fprintf( stdout, "\n" );
		for( k1=0; k1<lgth1; k1++ )
		{
			fprintf( stdout, "%d", k1 );
			for( k2=0; k2<lgth2; k2++ )
				fprintf( stdout, "%2.1f ", impmtx[k1][k2] );
			fprintf( stdout, "\n" );
		}
		exit( 1 );
	}
#endif
}
