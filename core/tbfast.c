//done
#include "mltaln.h"

#define DEBUG 0
#define IODEBUG 0
#define SCOREOUT 0
#define SHISHAGONYU 0 // for debug

static int nadd;
static int treein;
static int topin;
static int treeout;
static int distout;
static int noalign;
static int multidist;
static int subalignment;
static int subalignmentoffset;
static int keeplength;
static int ndeleted;
static int mapout;
static int smoothing;
static int specifictarget;
static int callpairlocalalign;
static int outputhat23;
static int compacttree = 0;

typedef struct _jobtable
{
    int i;  
    int j;  
} Jobtable;

typedef struct _msacompactdistmtxthread_arg // single thread demo tsukau
{
	int njob;
	int thread_no;
	int *selfscore;
	double **partmtx;
	char **seq;
	int **skiptable;
	double *mindist;
	int *mindistfrom;
 	int *jobpospt;
#ifdef enablemultithread
	pthread_mutex_t *mutex;
#endif
} msacompactdistmtxthread_arg_t;

#ifdef enablemultithread
typedef struct _distancematrixthread_arg
{
	int njob;
	int thread_no;
	int *selfscore;
	double **iscore;
	char **seq;
	int **skiptable;
	Jobtable *jobpospt;
	pthread_mutex_t *mutex;
} distancematrixthread_arg_t;

typedef struct _treebasethread_arg
{
	int thread_no;
	int *nrunpt;
	int njob;
	int *nlen;
	int *jobpospt;
	int ***topol;
	Treedep *dep;
	char **aseq;
	double *effarr;
	int *alloclenpt;
	LocalHom **localhomtable;
	RNApair ***singlerna;
	double *effarr_kozo;
	int *fftlog;
	char *mergeoralign;
	int *targetmap;
	pthread_mutex_t *mutex;
	pthread_cond_t *treecond;
} treebasethread_arg_t;
#endif

static void arguments( int argc, char *argv[], int *pac, char **pav, int *tac, char **tav ) // 2 kai yobaremasu.
{
    int c;
	int i;

	nthread = 1;
	outnumber = 0;
	scoreout = 0;
	spscoreout = 0;
	treein = 0;
	topin = 0;
	rnaprediction = 'm';
	rnakozo = 0;
	nevermemsave = 0;
	inputfile = NULL;
	addfile = NULL;
	addprofile = 1;
	fftkeika = 0;
	constraint = 0;
	nblosum = 62;
	fmodel = 0;
	calledByXced = 0;
	devide = 0;
	use_fft = 0; // chuui
	force_fft = 0;
	fftscore = 1;
	fftRepeatStop = 0;
	fftNoAnchStop = 0;
    weight = 3;
    utree = 1;
	tbutree = 1;
    refine = 0;
    check = 1;
    cut = 0.0;
    disp = 0;
    outgap = 1;
    alg = 'A';
    mix = 0;
	tbitr = 0;
	scmtd = 5;
	tbweight = 0;
	tbrweight = 3;
	checkC = 0;
	treemethod = 'X';
	sueff_global = 0.1;
	contin = 0;
	scoremtx = 1;
	kobetsubunkatsu = 0;
//	dorp = NOTSPECIFIED;
	ppenalty_dist = NOTSPECIFIED;
	ppenalty = NOTSPECIFIED;
	penalty_shift_factor = 1000.0;
	ppenalty_ex = NOTSPECIFIED;
	poffset = NOTSPECIFIED;
	kimuraR = NOTSPECIFIED;
	pamN = NOTSPECIFIED;
	geta2 = GETA2;
	fftWinSize = NOTSPECIFIED;
	fftThreshold = NOTSPECIFIED;
	RNAppenalty = NOTSPECIFIED;
	RNAppenalty_ex = NOTSPECIFIED;
	RNApthr = NOTSPECIFIED;
	TMorJTT = JTT;
	consweight_multi = 1.0;
	consweight_rna = 0.0;
	multidist = 0;
	subalignment = 0;
	subalignmentoffset = 0;
	legacygapcost = 0;
	specificityconsideration = 0.0;
	keeplength = 0;
	mapout = 0;
	smoothing = 0;
	specifictarget = 0;
	callpairlocalalign = 0;
	outputhat23 = 0;
	nwildcard = 0;

	if( pac )
	{
		pav[0] = "tbfast-pair";
		*pac = 1;
		tav[0] = "tbfast";
		*tac = 1;
	
		for( i=0; i<argc; i++ )
		{
			if( argv[i][0] == '_' )
			{
				callpairlocalalign = 1;
//				reporterr( "start\n" );
	
				for( i++; i<argc; i++ )
				{
					if( argv[i][0] == '_' )
					{
//						reporterr( "end\n" );
						goto pavend;
					}
					pav[*pac] = argv[i];
					*pac += 1;
//					reporterr( "%s\n", argv[i] );
				}
			}
		}
	
	
		pavend:

//		reporterr( "i=%d\n", i );
		for( i++; i<argc; i++ )
		{
			tav[*tac] = argv[i];
			*tac += 1;
		}

		argc -= *pac + 1;
		argv += *pac + 1;
	
//		reporterr( "argc in tbfast = %d\n", argc );
//		reporterr( "*pac in tbfast = %d\n", *pac );
//		for( i=0; i<*tac; i++ ) reporterr( "%s\n", tav[i] );
	}
	else
	{
//		reporterr( "SECOND TIME\n" );
	}

//	reporterr( "*argv = %s\n", *argv );

    while( --argc > 0 && (*++argv)[0] == '-' )
	{
//		reporterr( "(*argv)[0] = %s\n", (*argv) );
        while ( ( c = *++argv[0] ) )
		{
//			reporterr( "c=%c\n", c );
            switch( c )
            {
				case 'i':
					inputfile = *++argv;
//					fprintf( stderr, "inputfile = %s\n", inputfile );
					--argc;
                    goto nextoption;
				case 'I':
					nadd = myatoi( *++argv );
//					fprintf( stderr, "nadd = %d\n", nadd );
					--argc;
					goto nextoption;
				case 'e':
					RNApthr = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'o':
					RNAppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
					--argc;
					goto nextoption;
				case 'V':
					ppenalty_dist = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'f':
					ppenalty = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty = %d\n", ppenalty );
					--argc;
					goto nextoption;
				case 'Q':
					penalty_shift_factor = atof( *++argv );
					--argc;
					goto nextoption;
				case 'g':
					ppenalty_ex = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "ppenalty_ex = %d\n", ppenalty_ex );
					--argc;
					goto nextoption;
				case 'h':
					poffset = (int)( atof( *++argv ) * 1000 - 0.5 );
//					fprintf( stderr, "poffset = %d\n", poffset );
					--argc;
					goto nextoption;
				case 'k':
					kimuraR = myatoi( *++argv );
//					fprintf( stderr, "kappa = %d\n", kimuraR );
					--argc;
					goto nextoption;
				case 'b':
					nblosum = myatoi( *++argv );
					scoremtx = 1;
//					fprintf( stderr, "blosum %d / kimura 200\n", nblosum );
					--argc;
					goto nextoption;
				case 'j':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = JTT;
//					fprintf( stderr, "jtt/kimura %d\n", pamN );
					--argc;
					goto nextoption;
				case 'm':
					pamN = myatoi( *++argv );
					scoremtx = 0;
					TMorJTT = TM;
//					fprintf( stderr, "tm %d\n", pamN );
					--argc;
					goto nextoption;
				case 'l':
					fastathreshold = atof( *++argv );
					constraint = 2;
					--argc;
					goto nextoption;
				case 'r':
					consweight_rna = atof( *++argv );
					rnakozo = 1;
					--argc;
					goto nextoption;
				case 'c':
					consweight_multi = atof( *++argv );
					--argc;
					goto nextoption;
				case 'C':
					nthread = myatoi( *++argv );
//					fprintf( stderr, "nthread = %d\n", nthread );
					--argc; 
#ifndef enablemultithread
					nthread = 0;
#endif
					goto nextoption;
				case 's':
					specificityconsideration = (double)myatof( *++argv );
//					fprintf( stderr, "specificityconsideration = %f\n", specificityconsideration );
					--argc; 
					goto nextoption;
				case 'R':
					rnaprediction = 'r';
#if 1
				case 'a':
					fmodel = 1;
					break;
#endif
				case 'K':
					addprofile = 0;
					break;
				case 'y':
					distout = 1;
					break;
				case 't':
					treeout = 1;
					break;
				case 'T':
					noalign = 1;
					break;
				case 'D':
					dorp = 'd';
					break;
				case 'P':
					dorp = 'p';
					break;
				case 'L':
					legacygapcost = 1;
					break;
#if 1
				case 'O':
					outgap = 0;
					break;
#else
				case 'O':
					fftNoAnchStop = 1;
					break;
#endif
#if 0
				case 'S' :
					scoreout = 1; // for checking parallel calculation
					break;
#else
				case 'S' :
					spscoreout = 1; // 2014/Dec/30, sp score
					break;
#endif
				case 'H':
					subalignment = 1;
					subalignmentoffset = myatoi( *++argv );
					--argc;
					goto nextoption;
#if 0
				case 'e':
					fftscore = 0;
					break;
				case 'r':
					fmodel = -1;
					break;
				case 'R':
					fftRepeatStop = 1;
					break;
				case 's':
					treemethod = 's';
					break;
#endif
				case 'X':
					treemethod = 'X';
					sueff_global = atof( *++argv );
//					fprintf( stderr, "sueff_global = %f\n", sueff_global );
					--argc;
					goto nextoption;
				case 'E':
					treemethod = 'E';
					break;
				case 'q':
					treemethod = 'q';
					break;
				case 'n' :
					outnumber = 1;
					break;
#if 0
				case 'a':
					alg = 'a';
					break;
				case 'H':
					alg = 'H';
					break;
				case 'Q':
					alg = 'Q';
					break;
#endif
				case '@':
					alg = 'd';
					break;
				case 'A':
					alg = 'A';
					break;
				case 'M':
					alg = 'M';
					break;
				case 'N':
					nevermemsave = 1;
					break;
				case 'B': // hitsuyou! memopt -M -B no tame
					break;
				case 'F':
					use_fft = 1;
					break;
				case 'G':
					force_fft = 1;
					use_fft = 1;
					break;
				case 'U':
					treein = 1;
					break;
#if 0
				case 'V':
					topin = 1;
					break;
#endif
				case 'u':
					tbrweight = 0;
					weight = 0;
					break;
				case 'v':
					tbrweight = 3;
					break;
				case 'd':
					multidist = 1;
					break;
#if 0
				case 'd':
					disp = 1;
					break;
#endif
/* Modified 01/08/27, default: user tree */
				case 'J':
					tbutree = 0;
					break;
/* modification end. */
				case 'z':
					fftThreshold = myatoi( *++argv );
					--argc; 
					goto nextoption;
				case 'w':
					fftWinSize = myatoi( *++argv );
					--argc;
					goto nextoption;
				case 'W':
					minimumweight = atof( *++argv );
//					fprintf( stderr, "minimumweight = %f\n", minimumweight );
					--argc;
					goto nextoption;
#if 0
				case 'Z':
					checkC = 1;
					break;
#endif
				case 'Y':
					keeplength = 1;
					break;
				case 'Z':
					mapout = 1;
					break;
				case 'p':
					smoothing = 1;
					break;
				case '=':
					specifictarget = 1;
					break;
				case ':':
					nwildcard = 1;
					break;
				case '+':
					outputhat23 = myatoi( *++argv );
//					fprintf( stderr, "outputhat23 = %f\n", outputhat23 );
					--argc;
					goto nextoption;
                default:
                    fprintf( stderr, "illegal option %c\n", c );
                    argc = 0;
                    break;
            }
		}
		nextoption:
			;
	}

//	reporterr( "argc=%d\n", argc );

    if( argc == 1 )
    {
        cut = atof( (*argv) );
        argc--;
    }
    if( argc != 0 ) 
    {
        fprintf( stderr, "tbfast options: Check source file !\n" );
        exit( 1 );
    }
	if( tbitr == 1 && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : o, m or u\n" );
		exit( 1 );
	}
	if( alg == 'C' && outgap == 0 )
	{
		fprintf( stderr, "conflicting options : C, o\n" );
		exit( 1 );
	}
}

#if 0
static void *distancematrixthread2( void *arg )
{
	distancematrixthread_arg_t *targ = (distancematrixthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	double *selfscore = targ->selfscore;
	double **iscore = targ->iscore;
	char **seq = targ->seq;
	Jobtable *jobpospt = targ->jobpospt;

	double ssi, ssj, bunbo;
	int i, j;

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		i = jobpospt->i;
		i++;
		if( i == njob-1 )
		{
			pthread_mutex_unlock( targ->mutex );
			return( NULL );
		}
		jobpospt->i = i;
		pthread_mutex_unlock( targ->mutex );

		ssi = selfscore[i];
		if( i % 10 == 0 ) fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
		for( j=i+1; j<njob; j++ )
		{
			ssj = selfscore[j];
			bunbo = MIN( ssi, ssj );
			if( bunbo == 0.0 )
				iscore[i][j-i] = 1.0;
			else
				iscore[i][j-i] = 1.0 - naivepairscore11( seq[i], seq[j], penalty_dist ) / bunbo;
		}
	}
}
#endif

static double preferenceval( int ori, int pos, int max ) // for debug
{
	pos -= ori;
	if( pos < 0 ) pos += max;
	return( 0.00000000000001 * pos );
}

static void *msacompactdisthalfmtxthread( void *arg ) // enablemultithread == 0 demo tsukau
{
	msacompactdistmtxthread_arg_t *targ = (msacompactdistmtxthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **partmtx = targ->partmtx;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	double *mindist = targ->mindist;
	int *mindistfrom = targ->mindistfrom;
 	int *jobpospt = targ->jobpospt;
	double tmpdist, preference, tmpdistx, tmpdisty;
	int i, j;

	while( 1 )
	{
#ifdef enablemultithread
		if( nthread ) 
		{
			pthread_mutex_lock( targ->mutex );
			i = *jobpospt;
			if( i == njob-1 )
			{
				pthread_mutex_unlock( targ->mutex );
				return( NULL );
			}
			*jobpospt = i+1;
			pthread_mutex_unlock( targ->mutex );
		}
		else
#endif
		{
			i = *jobpospt;
			if( i == njob-1 )
			{
				return( NULL );
			}
			*jobpospt = i+1;
		}

		if( i % 100 == 0 ) 
		{
			if( nthread )
				fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
			else
				fprintf( stderr, "\r% 5d / %d", i, njob );
		}

		for( j=i+1; j<njob; j++ ) 
		{
			tmpdist = distcompact_msa( seq[i], seq[j], skiptable[i], skiptable[j], selfscore[i], selfscore[j] ); // osoikedo,

			preference = preferenceval( i, j, njob );
			tmpdistx = tmpdist + preference;
			if( tmpdistx < mindist[i] )
			{
				mindist[i] = tmpdistx;
				mindistfrom[i] = j;
			}

			preference = preferenceval( j, i, njob );
			tmpdisty = tmpdist + preference;
			if( tmpdisty < mindist[j] )
			{
				mindist[j] = tmpdisty;
				mindistfrom[j] = i;
			}
			if( partmtx[i] ) partmtx[i][j] = tmpdist;
			if( partmtx[j] ) partmtx[j][i] = tmpdist;
		}
	}
}


#ifdef enablemultithread
#if 0
static void *distancematrixthread( void *arg ) // v7.2 ijou deha tsukawanaihazu
{
	distancematrixthread_arg_t *targ = (distancematrixthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	double *selfscore = targ->selfscore;
	double **iscore = targ->iscore;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	Jobtable *jobpospt = targ->jobpospt;

	double ssi, ssj, bunbo;
	int i, j;

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		j = jobpospt->j;
		i = jobpospt->i;
		j++;
		if( j == njob )
		{
			i++;
			j = i + 1;
			if( i == njob-1 )
			{
				pthread_mutex_unlock( targ->mutex );
				return( NULL );
			}
		}
		jobpospt->j = j;
		jobpospt->i = i;
		pthread_mutex_unlock( targ->mutex );


		if( j==i+1 && i % 10 == 0 ) fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );
		ssi = selfscore[i];
		ssj = selfscore[j];
		bunbo = MIN( ssi, ssj );
		if( bunbo == 0.0 )
			iscore[i][j-i] = 2.0; // 2013/Oct/17
		else
//			iscore[i][j-i] = ( 1.0 - naivepairscore11( seq[i], seq[j], penalty_dist ) / bunbo ) * 2.0; // 2013/Oct/17
			iscore[i][j-i] = ( 1.0 - naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) / bunbo ) * 2.0; // 2014/Aug/15 fast
		if( iscore[i][j-i] > 10 ) iscore[i][j-i] = 10.0; // 2015/Mar/17
	}
}
#else
static void *distancematrixthread( void *arg ) // v7.2 ijou deha tsukawanaihazu
{
	distancematrixthread_arg_t *targ = (distancematrixthread_arg_t *)arg;
	int njob = targ->njob;
	int thread_no = targ->thread_no;
	int *selfscore = targ->selfscore;
	double **iscore = targ->iscore;
	char **seq = targ->seq;
	int **skiptable = targ->skiptable;
	Jobtable *jobpospt = targ->jobpospt;

	int ssi, ssj, bunbo;
	int i, j;

	while( 1 )
	{
		pthread_mutex_lock( targ->mutex );
		i = jobpospt->i; // (jobpospt-i)++ dato, shuuryou hantei no mae ni ++ surunode, tomaranakunaru.


		if( i == njob-1 )
		{
			pthread_mutex_unlock( targ->mutex );
			return( NULL );
		}
		jobpospt->i += 1;
		pthread_mutex_unlock( targ->mutex );
		if( i % 100 == 0 ) fprintf( stderr, "\r% 5d / %d (thread %4d)", i, njob, thread_no );

		ssi = selfscore[i];
		for( j=i+1; j<njob; j++ )
		{
			ssj = selfscore[j];
			bunbo = MIN( ssi, ssj );
			if( bunbo == 0 )
				iscore[i][j-i] = 2.0; // 2013/Oct/17
			else
//			iscore[i][j-i] = ( 1.0 - naivepairscore11( seq[i], seq[j], penalty_dist ) / bunbo ) * 2.0; // 2013/Oct/17
				iscore[i][j-i] = ( 1.0 - naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) / bunbo ) * 2.0; // 2014/Aug/15 fast
			if( iscore[i][j-i] > 10.0 ) iscore[i][j-i] = 10.0; // 2015/Mar/17
		}
	}
}
#endif
static void *treebasethread( void *arg )
{
	treebasethread_arg_t *targ = (treebasethread_arg_t *)arg;
	int *nrunpt = targ->nrunpt;
	int thread_no = targ->thread_no;
	int njob = targ->njob;
	int *nlen = targ->nlen;
	int *jobpospt = targ->jobpospt;
	int ***topol = targ->topol;
	Treedep *dep = targ->dep;
	char **aseq = targ->aseq;
	double *effarr = targ->effarr;
	int *alloclen = targ->alloclenpt;
	LocalHom **localhomtable = targ->localhomtable;
	RNApair ***singlerna = targ->singlerna;
	double *effarr_kozo = targ->effarr_kozo;
	int *fftlog = targ->fftlog;
	int *targetmap = targ->targetmap;
	char *mergeoralign = targ->mergeoralign;

	char **mseq1, **mseq2;
	char **localcopy;
	int i, j, l;
	int len1, len2;
	int clus1, clus2;
	double pscore;
	char *indication1, *indication2;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	double *effarr1_kozo = NULL;
	double *effarr2_kozo = NULL;
	LocalHom ***localhomshrink = NULL;
	char *swaplist = NULL;
	int m1, m2;
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	RNApair ***grouprna1 = NULL, ***grouprna2 = NULL;
	double **dynamicmtx;


	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );
	localcopy = calloc( njob, sizeof( char * ) );
	dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );

	if( rnakozo && rnaprediction == 'm' )
	{
		grouprna1 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
		grouprna2 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
	}
	else
	{
		grouprna1 = grouprna2 = NULL;
	}

	effarr1 = AllocateDoubleVec( njob );
	effarr2 = AllocateDoubleVec( njob );
	indication1 = AllocateCharVec( 150 );
	indication2 = AllocateCharVec( 150 );
#if 0
#else
	swaplist = NULL;
	if( constraint )
	{
		if( specifictarget ) swaplist = calloc( njob, sizeof( char ) );
		localhomshrink = (LocalHom ***)calloc( njob, sizeof( LocalHom ** ) );
		for( i=0; i<njob; i++ )
			localhomshrink[i] = (LocalHom **)calloc( njob, sizeof( LocalHom *) );
	}
#endif
	effarr1_kozo = AllocateDoubleVec( njob ); //tsuneni allocate sareru.
	effarr2_kozo = AllocateDoubleVec( njob ); //tsuneni allocate sareru.
	for( i=0; i<njob; i++ ) effarr1_kozo[i] = 0.0;
	for( i=0; i<njob; i++ ) effarr2_kozo[i] = 0.0;


#if 0
#endif

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stderr, "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif

#if 0 //-> main thread
	if( constraint )
		calcimportance( njob, effarr, aseq, localhomtable );
#endif


//	writePre( njob, name, nlen, aseq, 0 );


//	for( l=0; l<njob-1; l++ )
	while( 1 )
	{

		pthread_mutex_lock( targ->mutex );
		l = *jobpospt;
		if( l == njob-1 )
		{
			pthread_mutex_unlock( targ->mutex );
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
			Falign_udpari_long( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL );
			A__align( NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1 );
			D__align( NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
			partA__align( NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL );
			G__align11( NULL, NULL, NULL, 0, 0, 0 ); // iru?
			free( mseq1 );
			free( mseq2 );
			free( localcopy );
			free( effarr1 );
			free( effarr2 );
			free( effarr1_kozo );
			free( effarr2_kozo );
			free( indication1 );
			free( indication2 );
			FreeDoubleMtx( dynamicmtx );
			if( rnakozo && rnaprediction == 'm' )
			{
				if( grouprna1 ) free( grouprna1 ); // nakami ha?
				if( grouprna2 ) free( grouprna2 ); // nakami ha?
				grouprna1 = grouprna2 = NULL;
			}
			if( constraint )
			{
				if( localhomshrink ) // nen no tame
				{
					for( i=0; i<njob; i++ )
					{
						free( localhomshrink[i] );
						localhomshrink[i] = NULL;
					}
					free( localhomshrink );
					localhomshrink = NULL;
				}
				if( specifictarget ) free( swaplist );
			}
			return( NULL );
		}
		*jobpospt = l+1;

		if( dep[l].child0 != -1 )
		{
			while( dep[dep[l].child0].done == 0 )
				pthread_cond_wait( targ->treecond, targ->mutex );
		}
		if( dep[l].child1 != -1 )
		{
			while( dep[dep[l].child1].done == 0 )
				pthread_cond_wait( targ->treecond, targ->mutex );
		}
//		while( *nrunpt >= nthread )
//			pthread_cond_wait( targ->treecond, targ->mutex );
		(*nrunpt)++;

//		pthread_mutex_unlock( targ->mutex );

		if( mergeoralign[l] == 'n' )
		{
//			fprintf( stderr, "SKIP!\n" );
			dep[l].done = 1;
			(*nrunpt)--;
			pthread_cond_broadcast( targ->treecond );
			free( topol[l][0] );
			free( topol[l][1] );
			free( topol[l] );
			pthread_mutex_unlock( targ->mutex );
			continue;
		}



		m1 = topol[l][0][0];
		m2 = topol[l][1][0];

//		fprintf( stderr, "\ndistfromtip = %f\n", dep[l].distfromtip );
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip - 0.5 );
		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );

//		pthread_mutex_lock( targ->mutex );



        len1 = strlen( aseq[m1] );
        len2 = strlen( aseq[m2] );
        if( *alloclen <= len1 + len2 )
        {
			fprintf( stderr, "\nReallocating (by thread %d) ..", thread_no );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  );
			fprintf( stderr, "done. *alloclen = %d\n", *alloclen );
		}

		for( i=0; (j=topol[l][0][i])!=-1; i++ )
		{
			localcopy[j] = calloc( *alloclen, sizeof( char ) );
			strcpy( localcopy[j], aseq[j] );
		}
		for( i=0; (j=topol[l][1][i])!=-1; i++ )
		{
			localcopy[j] = calloc( *alloclen, sizeof( char ) );
			strcpy( localcopy[j], aseq[j] );
		}

		pthread_mutex_unlock( targ->mutex );



		if( effarr_kozo )
		{
			clus1 = fastconjuction_noname_kozo( topol[l][0], localcopy, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1 );
			clus2 = fastconjuction_noname_kozo( topol[l][1], localcopy, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2 );
		}
#if 0
		else if( specifictarget )
		{
			clus1 = fastconjuction_target( topol[l][0], localcopy, mseq1, effarr1, effarr, indication1, minimumweight, targetmap );
			clus2 = fastconjuction_target( topol[l][1], localcopy, mseq2, effarr2, effarr, indication2, minimumweight, targetmap );
		}
#endif
		else
		{
			clus1 = fastconjuction_noname( topol[l][0], localcopy, mseq1, effarr1, effarr, indication1, minimumweight );
			clus2 = fastconjuction_noname( topol[l][1], localcopy, mseq2, effarr2, effarr, indication2, minimumweight );
		}

#if 1
		fprintf( stderr, "\rSTEP % 5d /%d (thread %4d) ", l+1, njob-1, thread_no );
#else
		fprintf( stderr, "STEP %d /%d (thread %d) \n", l+1, njob-1, thread_no );
		fprintf( stderr, "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, ", child1 = %d\n", dep[l].child0 );
		fprintf( stderr, "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, ", child2 = %d\n", dep[l].child1 );

		fprintf( stderr, "Group1's lengths = " );
		for( i=0; i<clus1; i++ ) fprintf( stderr, "%d ", strlen( mseq1[i] ) );
		fprintf( stderr, "\n" );
		fprintf( stderr, "Group2's lengths = " );
		for( i=0; i<clus2; i++ ) fprintf( stderr, "%d ", strlen( mseq2[i] ) );
		fprintf( stderr, "\n" );
		
#endif



//		for( i=0; i<clus1; i++ ) fprintf( stderr, "## STEP%d-eff for mseq1-%d %f\n", l+1, i, effarr1[i] );

		if( constraint )
		{
			if( specifictarget )
				fastshrinklocalhom_target( topol[l][0], topol[l][1], localhomtable, localhomshrink, swaplist, targetmap );
			else
				fastshrinklocalhom_half( topol[l][0], topol[l][1], localhomtable, localhomshrink );
//			msfastshrinklocalhom( topol[l][0], topol[l][1], localhomtable, localhomshrink );
//			fprintf( stderr, "localhomshrink =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
//			weightimportance4( clus1, clus2, effarr1, effarr2, localhomshrink );
//			fprintf( stderr, "after weight =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
		}
		if( rnakozo && rnaprediction == 'm' )
		{
			makegrouprna( grouprna1, singlerna, topol[l][0] );
			makegrouprna( grouprna2, singlerna, topol[l][1] );
		}


/*
		fprintf( stderr, "before align all\n" );
		display( localcopy, njob );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/

		if( !nevermemsave && ( constraint != 2  && alg != 'M'  && ( len1 > 30000 || len2 > 30000 ) ) )
		{
			fprintf( stderr, "\nlen1=%d, len2=%d, Switching to the memsave mode.\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}
		

//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000 );
		else						   ffttry = 0;
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000 ); // v6.708
//		fprintf( stderr, "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );
//		fprintf( stderr, "f=%d, clus1=%d, fftlog[m1]=%d, clus2=%d, fftlog[m2]=%d\n", ffttry, clus1, fftlog[m1], clus2, fftlog[m2] );
		if( constraint == 2 )
		{
			if( alg == 'M' )
			{
				fprintf( stderr, "\n\nMemory saving mode is not supported.\n\n" );
				exit( 1 );
			}
			fprintf( stderr, "c" );
			if( alg == 'A' )
			{
				imp_match_init_strict( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, swaplist, 1, topol[l][0], topol[l][1] );
				if( rnakozo ) imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, NULL, NULL, NULL );
				pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, -1, -1 );
			}
			if( alg == 'd' )
			{
				imp_match_init_strictD( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, swaplist, 1, topol[l][0], topol[l][1]  );
				if( rnakozo ) imp_rnaD( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, NULL, NULL, NULL );
				pscore = D__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
			}
			else if( alg == 'Q' )
			{
				fprintf( stderr, "Not supported\n" );
				exit( 1 );
			}
		}
		else if( force_fft || ( use_fft && ffttry ) )
		{
			fprintf( stderr, "f" );
			if( alg == 'M' )
			{
				fprintf( stderr, "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
		}
		else
		{
			fprintf( stderr, "d" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					fprintf( stderr, "m" );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					break;
				case( 'A' ):
					pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, -1, -1 );
					break;
				case( 'd' ):
					pscore = D__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}

		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

#if SCOREOUT
		fprintf( stderr, "score = %10.2f\n", pscore );
#endif

/*
		fprintf( stderr, "after align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "after align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/

//		writePre( njob, name, nlen, localcopy, 0 );

		if( disp ) display( localcopy, njob );




		pthread_mutex_lock( targ->mutex );
		dep[l].done = 1;
		(*nrunpt)--;
		pthread_cond_broadcast( targ->treecond );

		for( i=0; (j=topol[l][0][i])!=-1; i++ )
			strcpy( aseq[j], localcopy[j] );
		for( i=0; (j=topol[l][1][i])!=-1; i++ )
			strcpy( aseq[j], localcopy[j] );

		pthread_mutex_unlock( targ->mutex );

		for( i=0; (j=topol[l][0][i])!=-1; i++ )
			free( localcopy[j] );
		for( i=0; (j=topol[l][1][i])!=-1; i++ )
			free( localcopy[j] );
		free( topol[l][0] );
		free( topol[l][1] );
		free( topol[l] );

	}
}
#endif

void treebase( int *nlen, char **aseq, int nadd, char *mergeoralign, char **mseq1, char **mseq2, int ***topol, Treedep *dep, double *effarr, int *alloclen, LocalHom **localhomtable, RNApair ***singlerna, double *effarr_kozo, int *targetmap, int *targetmapr, int ntarget )
{
	int i, l, m;
	int len1nocommongap, len2nocommongap;
	int len1, len2;
	int clus1, clus2;
	double pscore, tscore;
	char *indication1, *indication2;
	double *effarr1 = NULL;
	double *effarr2 = NULL;
	double *effarr1_kozo = NULL;
	double *effarr2_kozo = NULL;
	LocalHom ***localhomshrink = NULL;
	char *swaplist = NULL;
	int *fftlog;
	int m1, m2;
	int *gaplen;
	int *gapmap;
	int *alreadyaligned;
//	double dumfl = 0.0;
	double dumdb = 0.0;
	int ffttry;
	RNApair ***grouprna1 = NULL, ***grouprna2 = NULL;
	static double **dynamicmtx;
	int gapmaplen;

	if( rnakozo && rnaprediction == 'm' )
	{
		grouprna1 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
		grouprna2 = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
	}
	else
	{
		grouprna1 = grouprna2 = NULL;
	}

	fftlog = AllocateIntVec( njob );
	effarr1 = AllocateDoubleVec( njob );
	effarr2 = AllocateDoubleVec( njob );
	indication1 = AllocateCharVec( 150 );
	indication2 = AllocateCharVec( 150 );
	gaplen = AllocateIntVec( *alloclen+10 );
	gapmap = AllocateIntVec( *alloclen+10 );
	alreadyaligned = AllocateIntVec( njob );
	dynamicmtx = AllocateDoubleMtx( nalphabets, nalphabets );
#if 0
#else
	swaplist = NULL;
	if( constraint )
	{
		if( specifictarget ) swaplist = calloc( njob, sizeof( char ) );
		localhomshrink = (LocalHom ***)calloc( njob, sizeof( LocalHom ** ) );
		for( i=0; i<njob; i++ )
			localhomshrink[i] = (LocalHom **)calloc( njob, sizeof( LocalHom *) );
	}
#endif
	effarr1_kozo = AllocateDoubleVec( njob ); //tsuneni allocate sareru.
	effarr2_kozo = AllocateDoubleVec( njob ); //tsuneni allocate sareru.
	for( i=0; i<njob; i++ ) effarr1_kozo[i] = 0.0;
	for( i=0; i<njob; i++ ) effarr2_kozo[i] = 0.0;

	for( i=0; i<njob-nadd; i++ ) alreadyaligned[i] = 1;
	for( i=njob-nadd; i<njob; i++ ) alreadyaligned[i] = 0;

	for( l=0; l<njob; l++ ) fftlog[l] = 1;

#if 0
	fprintf( stderr, "##### fftwinsize = %d, fftthreshold = %d\n", fftWinSize, fftThreshold );
#endif

#if 0
	for( i=0; i<njob; i++ )
		fprintf( stderr, "TBFAST effarr[%d] = %f\n", i, effarr[i] );
#endif



	if( constraint )
	{
//		calcimportance( njob, effarr, aseq, localhomtable );
//		dontcalcimportance( njob, effarr, aseq, localhomtable ); // CHUUIII!!!!!
		if( specifictarget )
			calcimportance_target( njob, ntarget, effarr, aseq, localhomtable, targetmap, targetmapr );
//			dontcalcimportance_target( njob, effarr, aseq, localhomtable, ntarget ); // CHUUIII!!!!!
		else
//			calcimportance( njob, effarr, aseq, localhomtable );
			calcimportance_half( njob, effarr, aseq, localhomtable );
	}


//	writePre( njob, name, nlen, aseq, 0 );

	tscore = 0.0;
	for( l=0; l<njob-1; l++ )
	{
//		fprintf( stderr, "\ndistfromtip = %f\n", dep[l].distfromtip );
		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, dep[l].distfromtip );
//		makedynamicmtx( dynamicmtx, n_dis_consweight_multi, ( dep[l].distfromtip - 0.2 ) * 3 );
		if( mergeoralign[l] == 'n' )
		{
//			fprintf( stderr, "SKIP!\n" );
			free( topol[l][0] );
			free( topol[l][1] );
			free( topol[l] );
			continue;
		}

		m1 = topol[l][0][0];
		m2 = topol[l][1][0];
        len1 = strlen( aseq[m1] );
        len2 = strlen( aseq[m2] );
        if( *alloclen < len1 + len2 )
        {
			fprintf( stderr, "\nReallocating.." );
			*alloclen = ( len1 + len2 ) + 1000;
			ReallocateCharMtx( aseq, njob, *alloclen + 10  );
			gaplen = realloc( gaplen, ( *alloclen + 10 ) * sizeof( int ) );
			if( gaplen == NULL )
			{
				fprintf( stderr, "Cannot realloc gaplen\n" );
				exit( 1 );
			}
			gapmap = realloc( gapmap, ( *alloclen + 10 ) * sizeof( int ) );
			if( gapmap == NULL )
			{
				fprintf( stderr, "Cannot realloc gapmap\n" );
				exit( 1 );
			}
			fprintf( stderr, "done. *alloclen = %d\n", *alloclen );
		}

		if( effarr_kozo )
		{
			clus1 = fastconjuction_noname_kozo( topol[l][0], aseq, mseq1, effarr1, effarr, effarr1_kozo, effarr_kozo, indication1 );
			clus2 = fastconjuction_noname_kozo( topol[l][1], aseq, mseq2, effarr2, effarr, effarr2_kozo, effarr_kozo, indication2 );
		}
#if 0
		else if( specifictarget )
		{
			clus1 = fastconjuction_target( topol[l][0], aseq, mseq1, effarr1, effarr, indication1, minimumweight, targetmap );
			clus2 = fastconjuction_target( topol[l][1], aseq, mseq2, effarr2, effarr, indication2, minimumweight, targetmap );
		}
#endif
		else
		{
			clus1 = fastconjuction_noname( topol[l][0], aseq, mseq1, effarr1, effarr, indication1, minimumweight );
			clus2 = fastconjuction_noname( topol[l][1], aseq, mseq2, effarr2, effarr, indication2, minimumweight );
		}

		if( mergeoralign[l] == '1' || mergeoralign[l] == '2' ) // only in serial version
		{
			newgapstr = "=";
		}
		else
			newgapstr = "-";


		len1nocommongap = len1;
		len2nocommongap = len2;
		if( mergeoralign[l] == '1' ) // nai
		{
			findcommongaps( clus2, mseq2, gapmap );
			commongappick( clus2, mseq2 );
			len2nocommongap = strlen( mseq2[0] );
		}
		else if( mergeoralign[l] == '2' )
		{
			findcommongaps( clus1, mseq1, gapmap );
			commongappick( clus1, mseq1 );
			len1nocommongap = strlen( mseq1[0] );
		}
		

		fprintf( trap_g, "\nSTEP-%d\n", l );
		fprintf( trap_g, "group1 = %s\n", indication1 );
		fprintf( trap_g, "group2 = %s\n", indication2 );

#if 1
		fprintf( stderr, "\rSTEP % 5d /%d ", l+1, njob-1 );
		fflush( stderr );
#else
		fprintf( stdout, "STEP %d /%d\n", l+1, njob-1 );
		fprintf( stderr, "STEP %d /%d\n", l+1, njob-1 );
		fprintf( stderr, "group1 = %.66s", indication1 );
		if( strlen( indication1 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, "\n" );
		fprintf( stderr, "group2 = %.66s", indication2 );
		if( strlen( indication2 ) > 66 ) fprintf( stderr, "..." );
		fprintf( stderr, "\n" );
#endif



//		for( i=0; i<clus1; i++ ) fprintf( stderr, "## STEP%d-eff for mseq1-%d %f\n", l+1, i, effarr1[i] );

		if( constraint )
		{
			if( specifictarget )
				fastshrinklocalhom_target( topol[l][0], topol[l][1], localhomtable, localhomshrink, swaplist, targetmap );
			else
				fastshrinklocalhom_half( topol[l][0], topol[l][1], localhomtable, localhomshrink );
//			msfastshrinklocalhom( topol[l][0], topol[l][1], localhomtable, localhomshrink );
//			fprintf( stderr, "localhomshrink =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
//			weightimportance4( clus1, clus2, effarr1, effarr2, localhomshrink );
//			fprintf( stderr, "after weight =\n" );
//			outlocalhompt( localhomshrink, clus1, clus2 );
		}
		if( rnakozo && rnaprediction == 'm' )
		{
			makegrouprna( grouprna1, singlerna, topol[l][0] );
			makegrouprna( grouprna2, singlerna, topol[l][1] );
		}



/*
		fprintf( stderr, "before align all\n" );
		display( aseq, njob );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "before align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/

		if( !nevermemsave && ( constraint != 2  && alg != 'M'  && ( len1 > 30000 || len2 > 30000 ) ) )
		{
			fprintf( stderr, "\nlen1=%d, len2=%d, Switching to the memsave mode.\n", len1, len2 );
			alg = 'M';
			if( commonIP ) FreeIntMtx( commonIP );
			commonIP = NULL;
			commonAlloc1 = 0;
			commonAlloc2 = 0;
		}
		

//		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 );
		if( fftlog[m1] && fftlog[m2] ) ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 1000 && clus2 < 1000 );
		else						   ffttry = 0;
//		ffttry = ( nlen[m1] > clus1 && nlen[m2] > clus2 && clus1 < 5000 && clus2 < 5000 ); // v6.708
//		fprintf( stderr, "f=%d, len1/fftlog[m1]=%f, clus1=%d, len2/fftlog[m2]=%f, clus2=%d\n", ffttry, (double)len1/fftlog[m1], clus1, (double)len2/fftlog[m2], clus2 );
//		fprintf( stderr, "f=%d, clus1=%d, fftlog[m1]=%d, clus2=%d, fftlog[m2]=%d\n", ffttry, clus1, fftlog[m1], clus2, fftlog[m2] );
		if( constraint == 2 )
		{
			if( alg == 'M' )
			{
				fprintf( stderr, "\n\nMemory saving mode is not supported.\n\n" );
				exit( 1 );
			}
			fprintf( stderr, "c" );
			if( alg == 'A' )
			{
				imp_match_init_strict( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, swaplist, 1, topol[l][0], topol[l][1]  );
				if( rnakozo ) imp_rna( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, NULL, NULL, NULL );
				pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, topol[l][0][0], 1 ); // reuse profiles
			}
			if( alg == 'd' )
			{
				imp_match_init_strictD( NULL, clus1, clus2, strlen( mseq1[0] ), strlen( mseq2[0] ), mseq1, mseq2, effarr1, effarr2, effarr1_kozo, effarr2_kozo, localhomshrink, swaplist, 1, topol[l][0], topol[l][1]  );
				if( rnakozo ) imp_rnaD( clus1, clus2, mseq1, mseq2, effarr1, effarr2, grouprna1, grouprna2, NULL, NULL, NULL );
				pscore = D__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, localhomshrink, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
			}
			else if( alg == 'Q' )
			{
				fprintf( stderr, "Not supported\n" );
				exit( 1 );
			}
		}
		else if( force_fft || ( use_fft && ffttry ) )
		{
			fprintf( stderr, "f" );
			if( alg == 'M' )
			{
				fprintf( stderr, "m" );
				pscore = Falign_udpari_long( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1 );
			}
			else
				pscore = Falign( NULL, NULL, dynamicmtx, mseq1, mseq2, effarr1, effarr2, NULL, NULL, clus1, clus2, *alloclen, fftlog+m1, NULL, 0, NULL );
		}
		else
		{
			fprintf( stderr, "d" );
			fftlog[m1] = 0;
			switch( alg )
			{
				case( 'a' ):
					pscore = Aalign( mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen );
					break;
				case( 'M' ):
					fprintf( stderr, "m" );
					pscore = MSalignmm( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					break;
				case( 'A' ):
					pscore = A__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap, topol[l][0][0], 1 ); // reuse profiles
					break;
				case( 'd' ):
					pscore = D__align( dynamicmtx, mseq1, mseq2, effarr1, effarr2, clus1, clus2, *alloclen, NULL, &dumdb, NULL, NULL, NULL, NULL, NULL, 0, NULL, outgap, outgap );
					break;
				default:
					ErrorExit( "ERROR IN SOURCE FILE" );
			}
		}

		nlen[m1] = 0.5 * ( nlen[m1] + nlen[m2] );

#if SCOREOUT
		fprintf( stderr, "score = %10.2f\n", pscore );
#endif
		tscore += pscore;
/*
		fprintf( stderr, "after align 1 %s \n", indication1 );
		display( mseq1, clus1 );
		fprintf( stderr, "\n" );
		fprintf( stderr, "after align 2 %s \n", indication2 );
		display( mseq2, clus2 );
		fprintf( stderr, "\n" );
*/

//		writePre( njob, name, nlen, aseq, 0 );

		if( disp ) display( aseq, njob );

		if( mergeoralign[l] == '1' ) // jissainiha nai. atarashii hairetsu ha saigo dakara.
		{
			reporterr( "Check source!!\n" );
			exit( 1 );
		}
		if( mergeoralign[l] == '2' )
		{
//			fprintf( stderr, ">mseq1[0] = \n%s\n", mseq1[0] );
//			fprintf( stderr, ">mseq2[0] = \n%s\n", mseq2[0] );
//			if( keeplength ) ndeleted += deletenewinsertions( clus1, clus2, mseq1, mseq2, NULL );
			gapmaplen = strlen( mseq1[0] )-len1nocommongap+len1;
			adjustgapmap( gapmaplen, gapmap, mseq1[0] );
			if( smoothing )
			{
				restorecommongapssmoothly( njob, njob-(clus1+clus2), aseq, topol[l][0], topol[l][1], gapmap, *alloclen, '-' );
				findnewgaps( clus1, 0, mseq1, gaplen );
				insertnewgaps_bothorders( njob, alreadyaligned, aseq, topol[l][0], topol[l][1], gaplen, gapmap, gapmaplen, *alloclen, alg, '-' );
			}
			else
			{
				restorecommongaps( njob, njob-(clus1+clus2), aseq, topol[l][0], topol[l][1], gapmap, *alloclen, '-' );
				findnewgaps( clus1, 0, mseq1, gaplen );
				insertnewgaps( njob, alreadyaligned, aseq, topol[l][0], topol[l][1], gaplen, gapmap, *alloclen, alg, '-' );
			}
#if 0
			for( i=0; i<njob; i++ ) eq2dash( aseq[i] );
			for( i=0; i<clus1; i++ ) eq2dash( mseq1[i] );
			for( i=0; i<clus2; i++ ) eq2dash( mseq2[i] );
#else
			eq2dashmatometehayaku( mseq1, clus1 );
			eq2dashmatometehayaku( mseq2, clus2 );
#endif
			for( i=0; (m=topol[l][1][i])>-1; i++ ) alreadyaligned[m] = 1;
		}

		free( topol[l][0] );
		free( topol[l][1] );
		free( topol[l] );
	}
	free( topol[l] );
#if SCOREOUT
	fprintf( stderr, "totalscore = %10.2f\n\n", tscore );
#endif
	if( rnakozo && rnaprediction == 'm' )
	{
		if( grouprna1 ) free( grouprna1 ); // nakami ha?
		if( grouprna2 ) free( grouprna2 ); // nakami ha?
		grouprna1 = grouprna2 = NULL;
	}
	if( constraint )
	{
		if( localhomshrink ) // nen no tame
		{
			for( i=0; i<njob; i++ )
			{
				free( localhomshrink[i] );
				localhomshrink[i] = NULL;
			}
			free( localhomshrink );
			localhomshrink = NULL;
		}
		if( specifictarget ) free( swaplist );
	}

	free( topol );
	free( fftlog );
	free( effarr1 );
	free( effarr2 );
	free( indication1 );
	free( indication2 );
	free( gaplen );
	free( gapmap );
	free( alreadyaligned );
	FreeDoubleMtx( dynamicmtx );
	free( effarr1_kozo );
	free( effarr2_kozo );
	Falign( NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, 0, NULL );
	D__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0 );
	A__align( NULL,  NULL, NULL, NULL, NULL, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, 0, 0, -1, -1 );
	imp_match_init_strictD( NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, NULL );
	imp_match_init_strict( NULL, 0, 0, 0, 0, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 0, NULL, NULL );
	FreeCommonIP();
}

static void WriteOptions( FILE *fp )
{

	if( dorp == 'd' ) fprintf( fp, "DNA\n" );
	else
	{
		if     ( scoremtx ==  0 ) fprintf( fp, "JTT %dPAM\n", pamN );
		else if( scoremtx ==  1 ) fprintf( fp, "BLOSUM %d\n", nblosum );
		else if( scoremtx ==  2 ) fprintf( fp, "M-Y\n" );
	}
    fprintf( stderr, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );
    if( use_fft ) fprintf( fp, "FFT on\n" );

	fprintf( fp, "tree-base method\n" );
	if( tbrweight == 0 ) fprintf( fp, "unweighted\n" );
	else if( tbrweight == 3 ) fprintf( fp, "clustalw-like weighting\n" );
	if( tbitr || tbweight ) 
	{
		fprintf( fp, "iterate at each step\n" );
		if( tbitr && tbrweight == 0 ) fprintf( fp, "  unweighted\n" ); 
		if( tbitr && tbrweight == 3 ) fprintf( fp, "  reversely weighted\n" ); 
		if( tbweight ) fprintf( fp, "  weighted\n" ); 
		fprintf( fp, "\n" );
	}

   	 fprintf( fp, "Gap Penalty = %+5.2f, %+5.2f, %+5.2f\n", (double)ppenalty/1000, (double)ppenalty_ex/1000, (double)poffset/1000 );

	if( alg == 'a' )
		fprintf( fp, "Algorithm A\n" );
	else if( alg == 'A' ) 
		fprintf( fp, "Algorithm A+\n" );
	else if( alg == 'C' ) 
		fprintf( fp, "Apgorithm A+/C\n" );
	else
		fprintf( fp, "Unknown algorithm\n" );

	if( treemethod == 'X' )
		fprintf( fp, "Tree = UPGMA (mix).\n" );
	else if( treemethod == 'E' )
		fprintf( fp, "Tree = UPGMA (average).\n" );
	else if( treemethod == 'q' )
		fprintf( fp, "Tree = Minimum linkage.\n" );
	else
		fprintf( fp, "Unknown tree.\n" );

    if( use_fft )
    {
        fprintf( fp, "FFT on\n" );
        if( dorp == 'd' )
            fprintf( fp, "Basis : 4 nucleotides\n" );
        else
        {
            if( fftscore )
                fprintf( fp, "Basis : Polarity and Volume\n" );
            else
                fprintf( fp, "Basis : 20 amino acids\n" );
        }
        fprintf( fp, "Threshold   of anchors = %d%%\n", fftThreshold );
        fprintf( fp, "window size of anchors = %dsites\n", fftWinSize );
    }
	else
        fprintf( fp, "FFT off\n" );
	fflush( fp );
}
	 
static double **preparepartmtx( int nseq )
{
	int i;
	double **val;
	double size;

	val = (double **)calloc( nseq, sizeof( double *) );;
	size = 0;

#if 0
	if( compacttree == 1 )
	{
		for( i=0; i<nseq; i++ )
		{
			size += (double)sizeof( double ) * nseq;
			if( size > maxdistmtxsize ) 
			{
				reporterr( "\n\nThe size of full distance matrix is estimated to exceed %.2fGB.\n", maxdistmtxsize / 1000 / 1000 /1000 );
				reporterr( "Will try the calculation using a %d x %d matrix.\n", nseq, i );
				reporterr( "This calculation will be slow due to the limited RAM space.\n", i, nseq );
				reporterr( "To avoid the slowdown, please try '--initialramusage xGB' (x>>%.2f),\n", maxdistmtxsize / 1000 / 1000 /1000 );
				reporterr( "if larger RAM space is available.\n" );
				reporterr( "Note that xGB is NOT the upper limit of RAM usage.\n" );
				reporterr( "Two to three times larger space may be used for building a guide tree.\n" );
				reporterr( "Memory usage of the MSA stage depends on similarity of input sequences.\n\n" );
//				reporterr( "If the RAM is small, try '--initialramusage xGB' with a smaller x value.\n" );
				reporterr( "The '--memsavetree' option uses smaller RAM space.\n" );
				reporterr( "If tree-like relationship can be ignored, try '--pileup' or '--randomchain'.\n\n" );
				reporterr( "The result of --initialramusage xGB is almost identical to the default, except for rounding differences.\n" );

				reporterr( "In the cases of --memsavetree, --pileup and --randomchain, the result differs from the default.\n\n" );
				break;
			}
			val[i] = (double *)calloc( nseq, sizeof( double ) );
		}
		if( i == nseq ) reporterr( "The full matrix will be used.\n" );

		for( ;i<nseq; i++ ) val[i] = NULL; // nen no tame
	}
	else
#endif
	{
		for( i=0; i<nseq; i++ ) val[i] = NULL; // nen no tame
	}
	return( val );
}


int main( int argc, char *argv[] )
{
	static int  *nlen = NULL;	
	static int *selfscore = NULL;
	int nogaplen;
	static char **name = NULL, **seq = NULL;
	static char **mseq1 = NULL, **mseq2 = NULL;
	static char **bseq = NULL;
	static double **iscore = NULL, **iscore_kozo = NULL;
	int **skiptable;
	static double *eff = NULL, *eff_kozo = NULL, *eff_kozo_mapped = NULL;
	int i, j, ien, ik, jk;
	static int ***topol = NULL, ***topol_kozo = NULL;
	static int *addmem;
	static Treedep *dep = NULL;
	static double **len = NULL, **len_kozo = NULL;
	FILE *prep = NULL;
	FILE *infp = NULL;
	FILE *orderfp = NULL;
	FILE *hat2p = NULL;
	double unweightedspscore;
	int alignmentlength;
	char *mergeoralign = NULL;
	int foundthebranch;
	int nsubalignments, maxmem;
	int **subtable;
	int *insubtable;
	int *preservegaps;
	char ***subalnpt;
	char *originalgaps = NULL;
	char **addbk = NULL;
	int **deletelist = NULL;
	FILE *dlf = NULL;
// for compacttree
	int *mindistfrom = NULL;
	double *mindist = NULL;
	double **partmtx = NULL;
// for compacttree

	char c;
	int alloclen;
	LocalHom **localhomtable = NULL;
	LocalHom *tmpptr;
	RNApair ***singlerna = NULL;
	double ssi, ssj, bunbo;
	static char *kozoarivec = NULL;
	int nkozo;
	int ntarget;
	int *targetmap = NULL, *targetmapr = NULL;
	int ilim, jst, jj;
	int pac, tac;
	char **pav, **tav;

	pav = calloc( argc, sizeof( char * ) );
	tav = calloc( argc, sizeof( char * ) );

	arguments( argc, argv, &pac, pav, &tac, tav );

	if( fastathreshold < 0.0001 ) constraint = 0; 

	if( inputfile )
	{
		infp = fopen( inputfile, "r" );
		if( !infp ) 
		{
			fprintf( stderr, "Cannot open %s\n", inputfile );
			exit( 1 );
		}
	}
	else    
		infp = stdin;

	getnumlen( infp );
	rewind( infp );



	nkozo = 0;

	if( njob < 2 )
	{
		fprintf( stderr, "At least 2 sequences should be input!\n"
						 "Only %d sequence found.\n", njob ); 
		exit( 1 );
	}

	if( subalignment )
	{
		readsubalignmentstable( njob, NULL, NULL, &nsubalignments, &maxmem );
		fprintf( stderr, "nsubalignments = %d\n", nsubalignments );
		fprintf( stderr, "maxmem = %d\n", maxmem );
		subtable = AllocateIntMtx( nsubalignments, maxmem+1 );
		insubtable = AllocateIntVec( njob );
		for( i=0; i<njob; i++ ) insubtable[i] = 0;
		preservegaps = AllocateIntVec( njob );
		for( i=0; i<njob; i++ ) preservegaps[i] = 0;
		subalnpt = AllocateCharCub( nsubalignments, maxmem, 0 );
		readsubalignmentstable( njob, subtable, preservegaps, NULL, NULL );
	}

	seq = AllocateCharMtx( njob, nlenmax+1 );
	mseq1 = AllocateCharMtx( njob, 0 );
	mseq2 = AllocateCharMtx( njob, 0 );

	name = AllocateCharMtx( njob, B+1 );
	nlen = AllocateIntVec( njob );
	selfscore = AllocateIntVec( njob );

	topol = AllocateIntCub( njob, 2, 0 );
	len = AllocateFloatMtx( njob, 2 );
	eff = AllocateDoubleVec( njob );
	kozoarivec = AllocateCharVec( njob );

	mergeoralign = AllocateCharVec( njob );

	dep = (Treedep *)calloc( njob, sizeof( Treedep ) );
	if( nadd ) addmem = AllocateIntVec( nadd+1 );


	if( tbutree ) iscore = AllocateFloatHalfMtx( njob ); // tbutree=0 no toki aln kara mtx wo keisan, compacttree dehanaitoki nomi iscore shiyou.

	ndeleted = 0;

#if 0
	readData( infp, name, nlen, seq );
#else
	readData_pointer( infp, name, nlen, seq );
	fclose( infp );
#endif

	if( specifictarget )
	{
		targetmap = calloc( njob, sizeof( int ) );
		ntarget = 0;
		for( i=0; i<njob; i++ )
		{
			targetmap[i] = -1;
			if( !strncmp( name[i]+1, "_focus_", 7 ) )
				targetmap[i] = ntarget++;
		}
		targetmapr = calloc( ntarget, sizeof( int ) );
		for( i=0; i<njob; i++ )
			if( targetmap[i] != -1 ) targetmapr[targetmap[i]] = i;

	}
	else
	{
		ntarget = njob;
		targetmap = calloc( njob, sizeof( int ) );
		targetmapr = calloc( njob, sizeof( int ) );
		for( i=0; i<njob; i++ )
			targetmap[i] = targetmapr[i] = i;
	}

#if 0
	for( i=0; i<njob; i++ )
		reporterr( "targetmap[%d] = %d\n", i, targetmap[i] );
	for( i=0; i<ntarget; i++ )
		reporterr( "targetmapr[%d] = %d\n", i, targetmapr[i] );
#endif

//	if( constraint && !noalign ) // 2016mar15 noalign tsuika
	if( constraint ) // 2016Jul31 noalign no toki no shori (l=0.0) ha mafft.tmpl ni idou
	{

		ilim = njob;
		localhomtable = (LocalHom **)calloc( ntarget, sizeof( LocalHom *) );
		for( i=0; i<ntarget; i++ )
		{
			localhomtable[i] = (LocalHom *)calloc( ilim, sizeof( LocalHom ) );
			for( j=0; j<ilim; j++ )
			{
				localhomtable[i][j].start1 = -1;
				localhomtable[i][j].end1 = -1;
				localhomtable[i][j].start2 = -1;
				localhomtable[i][j].end2 = -1;
				localhomtable[i][j].overlapaa = -1.0;
				localhomtable[i][j].opt = -1.0;
				localhomtable[i][j].importance = -1.0;
				localhomtable[i][j].next = NULL;
				localhomtable[i][j].nokori = 0;
				localhomtable[i][j].extended = -1;
				localhomtable[i][j].last = localhomtable[i]+j;
				localhomtable[i][j].korh = 'h';
			}
			if( !specifictarget ) ilim--;
		}

//		reporterr( "pac=%d\n", pac );
//		reporterr( "pav[0]=%s\n", pav[0] );
		if( callpairlocalalign )
		{
			pairlocalalign( njob, nlenmax, name, seq, iscore, localhomtable, pac, pav );
			arguments( tac, tav, NULL, NULL, NULL, NULL ); // anzen no tame
			callpairlocalalign = 1; // wakarinikui.
			if( fastathreshold < 0.0001 ) constraint = 0; 
//			fprintf( stderr, "blosum %d / kimura 200\n", nblosum );
//			fprintf( stderr, "scoremtx=%d\n", scoremtx );
//			fprintf( stderr, "fastathreshold=%f\n", fastathreshold );
//			fprintf( stderr, "constraing=%d\n", constraint );
//exit( 1 );
			for( ilim=njob, i=0; i<ntarget; i++ )
			{
				for( j=0; j<ilim; j++ )
				{
					for( tmpptr=localhomtable[i]+j; tmpptr; tmpptr=tmpptr->next )
					{
						if( tmpptr->opt == -1.0 ) continue;
#if SHISHAGONYU // for debug
						char buff[100];
						sprintf( buff, "%10.5f", tmpptr->opt );
						tmpptr->opt = 0.0;
						sscanf( buff, "%lf", &(tmpptr->opt) );
#endif
						tmpptr->opt = ( tmpptr->opt ) / 5.8  * 600;
					}
				}
				if( !specifictarget ) ilim--;
			}

			prep = fopen( "hat3.seed", "r" );
			if( prep )
			{
				fprintf( stderr, "Loading 'hat3.seed' ... " );
				if( specifictarget ) readlocalhomtable2_target( prep, njob, localhomtable, kozoarivec, targetmap ); // uwagakisarerukara koredehadame.
				else readlocalhomtable2_half( prep, njob, localhomtable, kozoarivec );                                       // uwagakisarerukara koredehadame.
				fclose( prep );
				fprintf( stderr, "\ndone.\n" );
			}
			else
				fprintf( stderr, "No hat3.seed.\n" );

			if( outputhat23 )
			{
				prep = fopen( "hat3", "w" );
				if( !prep ) ErrorExit( "Cannot open hat3 to write." );

				fprintf( stderr, "Writing hat3 for iterative refinement\n" );
				if( specifictarget )
					ilim = ntarget;
				else
					ilim = njob-1;	
				for( i=0; i<ilim; i++ ) 
				{
					if( specifictarget )
					{
						jst = 0;
						jj = 0;
					}
					else
					{
						jst = i;
						jj = 0;
					}
					for( j=jst; j<njob; j++, jj++ )
					{
						for( tmpptr=localhomtable[i]+jj; tmpptr; tmpptr=tmpptr->next )
						{
							if( tmpptr->opt == -1.0 ) continue;
							if( targetmap[j] == -1 || targetmap[i] < targetmap[j] )
								fprintf( prep, "%d %d %d %7.5f %d %d %d %d %c\n", targetmapr[i], j, tmpptr->overlapaa, tmpptr->opt/600*5.8, tmpptr->start1, tmpptr->end1, tmpptr->start2, tmpptr->end2, tmpptr->korh );
						}
					}
				}
				fclose( prep );

				prep = fopen( "hat2", "w" );
				WriteFloatHat2_pointer_halfmtx( prep, njob, name, iscore );
				fclose( prep );
			}
			else if( distout ) // choufuku shiterukedo, muda deha nai.
			{
				prep = fopen( "hat2", "w" );
				WriteFloatHat2_pointer_halfmtx( prep, njob, name, iscore );
				fclose( prep );
			}
		}
		else
		{
			fprintf( stderr, "Loading 'hat3' ... " );
			prep = fopen( "hat3", "r" );
			if( prep == NULL ) ErrorExit( "Make hat3." );
			if( specifictarget ) readlocalhomtable2_target( prep, njob, localhomtable, kozoarivec, targetmap );
			else readlocalhomtable2_half( prep, njob, localhomtable, kozoarivec );
			fclose( prep );
			fprintf( stderr, "\ndone.\n" );
		}



		nkozo = 0;
		for( i=0; i<njob; i++ ) 
		{
//			fprintf( stderr, "kozoarivec[%d] = %d\n", i, kozoarivec[i] );
			if( kozoarivec[i] ) nkozo++;
		}
		if( nkozo )
		{
			topol_kozo = AllocateIntCub( nkozo, 2, 0 );
			len_kozo = AllocateFloatMtx( nkozo, 2 );
			iscore_kozo = AllocateFloatHalfMtx( nkozo );
			eff_kozo = AllocateDoubleVec( nkozo );
			eff_kozo_mapped = AllocateDoubleVec( njob );
		}


#if 0
		if( specifictarget )
			outlocalhom_target( localhomtable, ntarget, njob );
		else
			outlocalhom_half( localhomtable, njob );
		exit( 1 );
#endif

#if 0
		fprintf( stderr, "Extending localhom ... " );
		extendlocalhom2( njob, localhomtable );
		fprintf( stderr, "done.\n" );
#endif
	}
	else
	{

		if( callpairlocalalign )
		{
			pairlocalalign( njob, nlenmax, name, seq, iscore, NULL, pac, pav );
			arguments( tac, tav, NULL, NULL, NULL, NULL ); // anzen no tame
			callpairlocalalign = 1; // wakarinikui.
			if( fastathreshold < 0.0001 ) constraint = 0; 
			fprintf( stderr, "blosum %d / kimura 200\n", nblosum );
			fprintf( stderr, "scoremtx=%d\n", scoremtx );
			fprintf( stderr, "fastathreshold=%f\n", fastathreshold );
		}
		if( distout || outputhat23 )
		{
			reporterr( "\nwriting hat2 (1)\n" );
			prep = fopen( "hat2", "w" );
			WriteFloatHat2_pointer_halfmtx( prep, njob, name, iscore );
			fclose( prep );
		}
	}


	free( tav );
	free( pav );
	constants( njob, seq );


#if 0
	fprintf( stderr, "params = %d, %d, %d\n", penalty, penalty_ex, offset );
#endif

	initSignalSM();

	initFiles();

	WriteOptions( trap_g );

	if( distout && !treeout && noalign )  // 2016Jul31. Free ha mada fukanzen.
	{
		writeData_pointer( prep_g, njob, name, nlen, seq );
		fprintf( stderr, "\n" ); 
		SHOWVERSION;
		goto chudan;
	}




	c = seqcheck( seq );
	if( c )
	{
		fprintf( stderr, "Illegal character %c\n", c );
		exit( 1 );
	}

//	writePre( njob, name, nlen, seq, 0 );

	if( treein )
	{
		int dumx, dumy;
		double dumz;
		treein = check_guidetreefile( &dumx, &dumy, &dumz );
		if( treein == 'C' )
		{
			compacttree = 2;
			treein = 0;
			use_fft = 0; // kankeinai?
		}

		// else treein = 1 no mama
	}

	reporterr( "treein = %d\n", treein );
	reporterr( "compacttree = %d\n", compacttree );

	if( nadd && keeplength )
	{
		originalgaps = (char *)calloc( nlenmax+1, sizeof( char) );
		recordoriginalgaps( originalgaps, njob-nadd, seq );

		if( mapout )
		{
			addbk = (char **)calloc( nadd+1, sizeof( char * ) );
			for( i=0; i<nadd; i++ )
			{
				ien = strlen( seq[njob-nadd+i] );
				addbk[i] = (char *)calloc( ien + 1, sizeof( char ) );
				gappick0( addbk[i], seq[njob-nadd+i] );
			}
			addbk[nadd] = NULL;
		}
		else
			addbk = NULL;
	}
	else
	{
		originalgaps = NULL;
		addbk = NULL;
	}

	if( treein )
	{
#if 0
		if( nkozo )
		{
			fprintf( stderr, "Both structure and user tree have been given. Not yet supported!\n" );
			exit( 1 );
		}
#endif
		loadtree( njob, topol, len, name, nlen, dep, treeout );
//		loadtop( njob, topol, len, name, NULL, dep ); // 2015/Jan/13, not yet checked
		fprintf( stderr, "\ndone.\n\n" );
	}
	else
	{
		if( tbutree == 0 && compacttree ) // compacttree no toki ha treein ha 0 de uwagaki sarete iru.
		{
			iscore = NULL;// tsukawanai
			reporterr( "Making a compact tree from msa, step 1.. \n" );
			skiptable = AllocateIntMtx( njob, 0 );
			makeskiptable( njob, skiptable, seq ); // allocate suru.
			mindistfrom = (int *)calloc( njob, sizeof( int ) );
			mindist = (double *)calloc( njob, sizeof( double ) );
			partmtx = preparepartmtx( njob );

			for( i=0; i<njob; i++ ) // disttbfast deha kokoniha nakatta.
			{
//				selfscore[i] = naivepairscore11( seq[i], seq[i], penalty_dist );
				selfscore[i] = (int)naivepairscorefast( seq[i], seq[i], skiptable[i], skiptable[i], penalty_dist );
//				fprintf( stderr, "penalty = %d\n", penalty );
//				fprintf( stderr, "penalty_dist = %d\n", penalty_dist );
			}
#ifdef enablemultithread
			if( nthread > 0 )
			{
				msacompactdistmtxthread_arg_t *targ;
				int jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;
				double **mindistthread;
				int **mindistfromthread;

				mindistthread = AllocateDoubleMtx( nthread, njob );
				mindistfromthread = AllocateIntMtx( nthread, njob );
				targ = calloc( nthread, sizeof( msacompactdistmtxthread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex, NULL );
				jobpos = 0;

				for( i=0; i<nthread; i++ )
				{
					for( j=0; j<njob; j++ )
					{
						mindistthread[i][j] = 999.9;
						mindistfromthread[i][j] = -1;
					}
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].selfscore = selfscore;
					targ[i].partmtx = partmtx;
					targ[i].seq = seq;
					targ[i].skiptable = skiptable;
					targ[i].jobpospt = &jobpos;
					targ[i].mindistfrom = mindistfromthread[i];
					targ[i].mindist = mindistthread[i];
					targ[i].mutex = &mutex;
	
					pthread_create( handle+i, NULL, msacompactdisthalfmtxthread, (void *)(targ+i) );
				}
	
				for( i=0; i<nthread; i++ ) pthread_join( handle[i], NULL );
				pthread_mutex_destroy( &mutex );

				for( i=0; i<njob; i++ )
				{
					mindist[i] = 999.9;
					mindistfrom[i] = -1;
					for( j=0; j<nthread; j++ )
					{
						if( mindistthread[j][i] < mindist[i] )
						{
							mindist[i] = mindistthread[j][i];
							mindistfrom[i] = mindistfromthread[j][i];
						}
					}
				}
				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug

				free( handle );
				free( targ );
				FreeDoubleMtx( mindistthread );
				FreeIntMtx( mindistfromthread );
			}
			else
#endif
			{
				msacompactdistmtxthread_arg_t *targ;
				int jobpos;
				jobpos = 0;
				targ = calloc( 1, sizeof( msacompactdistmtxthread_arg_t ) );

				{
					for( j=0; j<njob; j++ )
					{
						mindist[j] = 999.9;
						mindistfrom[j] = -1;
					}
					targ[0].thread_no = 0;
					targ[0].njob = njob;
					targ[0].selfscore = selfscore;
					targ[0].partmtx = partmtx;
					targ[0].seq = seq;
					targ[0].skiptable = skiptable;
					targ[0].jobpospt = &jobpos;
					targ[0].mindistfrom = mindistfrom;
					targ[0].mindist = mindist;
	
					msacompactdisthalfmtxthread( targ );
//					msacompactdistmtxthread( targ );
				}
				free( targ );
				for( i=0; i<njob; i++ ) mindist[i] -= preferenceval( i, mindistfrom[i], njob ); // for debug
			}
//			free( selfscore ); selfscore = NULL; // mada tsukau
//			FreeCharMtx( bseq ); bseq = NULL; // mada tsukau
//			if( skiptable) FreeIntMtx( skiptable ); skiptable = NULL;

//			for( i=0; i<njob; i++ ) printf( "mindist[%d] = %f\n", i, mindist[i] );
//			exit( 1 );
			reporterr( "\rdone.                                          \n" );
		}
		else if( tbutree == 0 && compacttree == 0 )
		{
			reporterr( "Making a distance matrix from msa .. \n" );
//			reporterr( "Bug.  This function should not be used in versions >=7.2.  Please email kazutaka.katoh@aist.go.jp\n" );
//			fflush( stderr );
//			exit( 1 );
			iscore = AllocateFloatHalfMtx( njob );  // tbutree == 0 no baai ha allocate sareteinainode

			for( i=1; i<njob; i++ ) 
			{
				if( nlen[i] != nlen[0] ) 
				{
					fprintf( stderr, "Input pre-aligned seqences or make hat2.\n" );
					exit( 1 );
				}
			}
	
			skiptable = AllocateIntMtx( njob, 0 );
			makeskiptable( njob, skiptable, seq ); // allocate suru.
			ien = njob-1;
			for( i=0; i<njob; i++ ) 
			{
//				selfscore[i] = naivepairscore11( seq[i], seq[i], penalty_dist );
				selfscore[i] = (int)naivepairscorefast( seq[i], seq[i], skiptable[i], skiptable[i], penalty_dist );
//				fprintf( stderr, "penalty = %d\n", penalty );
//				fprintf( stderr, "penalty_dist = %d\n", penalty_dist );
			}
#ifdef enablemultithread
			if( nthread > 0 )
			{
				distancematrixthread_arg_t *targ;
				Jobtable jobpos;
				pthread_t *handle;
				pthread_mutex_t mutex;

				jobpos.i = 0;
				jobpos.j = 0;

				targ = calloc( nthread, sizeof( distancematrixthread_arg_t ) );
				handle = calloc( nthread, sizeof( pthread_t ) );
				pthread_mutex_init( &mutex, NULL );

				for( i=0; i<nthread; i++ )
				{
					targ[i].thread_no = i;
					targ[i].njob = njob;
					targ[i].selfscore = selfscore;
					targ[i].iscore = iscore;
					targ[i].seq = seq;
					targ[i].skiptable = skiptable;
					targ[i].jobpospt = &jobpos;
					targ[i].mutex = &mutex;

					pthread_create( handle+i, NULL, distancematrixthread, (void *)(targ+i) );
				}

				for( i=0; i<nthread; i++ )
				{
					pthread_join( handle[i], NULL );
				}
				pthread_mutex_destroy( &mutex );
				free( handle );
				free( targ );
			}
			else
#endif
			{
				for( i=0; i<ien; i++ ) 
				{
					if( i % 10 == 0 )
					{
						fprintf( stderr, "\r% 5d / %d", i, ien );
						fflush( stderr );
					}
					ssi = selfscore[i];
					for( j=i+1; j<njob; j++ ) 
					{
						ssj = selfscore[j];
						bunbo = MIN( ssi, ssj );
						if( bunbo == 0.0 )
							iscore[i][j-i] = 2.0; // 2013/Oct/17 2bai
						else
//							iscore[i][j-i] = 1.0 - naivepairscore11( seq[i], seq[j], penalty_dist ) / MIN( selfscore[i], selfscore[j] );
//							iscore[i][j-i] = ( 1.0 - naivepairscore11( seq[i], seq[j], penalty_dist ) / bunbo ) * 2.0; // 2013/Oct/17 2bai
							iscore[i][j-i] = ( 1.0 - naivepairscorefast( seq[i], seq[j], skiptable[i], skiptable[j], penalty_dist ) / bunbo ) * 2.0; // 2014/Aug/15 fast
						if( iscore[i][j-i] > 10 ) iscore[i][j-i] = 10.0; // 2015/Mar/17
//exit( 1 );
		
#if 0
						fprintf( stderr, "### ssj = %f\n", ssj );
						fprintf( stderr, "### selfscore[i] = %f\n", selfscore[i] );
						fprintf( stderr, "### selfscore[j] = %f\n", selfscore[j] );
						fprintf( stderr, "### rawscore = %f\n", naivepairscore11( seq[i], seq[j], penalty_dist ) );
#endif
					}
				}
			}
//			fprintf( stderr, "\ndone.\n\n" );
			FreeIntMtx( skiptable );
//			fflush( stderr );
			reporterr( "\rdone.                                           \n" );

		}
		else
		{
			if( callpairlocalalign )
			{
				if( multidist )
				{
					reporterr( "Bug in v7.290.  Please email kazutaka.katoh@aist.go.jp\n" );
					exit( 1 );
				}
#if 0
				prep = fopen( "hat2", "w" );
				if( !prep ) ErrorExit( "Cannot open hat2." );
				WriteFloatHat2_pointer_halfmtx( prep, njob, name, iscore ); // jissiha double
				fclose( prep );
#endif
			}
			else
			{
				if( multidist )
				{
					fprintf( stderr, "Loading 'hat2n' (aligned sequences - new sequences) ... " );
					prep = fopen( "hat2n", "r" );
					if( prep == NULL ) ErrorExit( "Make hat2." );
					readhat2_doublehalf_pointer( prep, njob, name, iscore );
					fclose( prep );
					fprintf( stderr, "done.\n" );
				
					fprintf( stderr, "Loading 'hat2i' (aligned sequences) ... " );
					prep = fopen( "hat2i", "r" );
					if( prep == NULL ) ErrorExit( "Make hat2i." );
					readhat2_doublehalf_pointer( prep, njob-nadd, name, iscore );
					fclose( prep );
					fprintf( stderr, "done.\n" );
				}
				else
				{
					fprintf( stderr, "Loading 'hat2' ... " );
					prep = fopen( "hat2", "r" );
					if( prep == NULL ) ErrorExit( "Make hat2." );
					readhat2_doublehalf_pointer( prep, njob, name, iscore );
					fclose( prep );
					fprintf( stderr, "done.\n" );
				}

				if( distout ) // callpairlocalalign == 1 no toki ha ue de shorizumi.
				{
					reporterr( "\nwriting hat2 (2)\n" );
					hat2p = fopen( "hat2", "w" );
					WriteFloatHat2_pointer_halfmtx( hat2p, njob, name, iscore );
					fclose( hat2p );
				}
			}
//			for( i=0; i<njob-1; i++ ) for( j=i+1; j<njob; j++ ) printf( "dist %d-%d = %f\n", i, j, iscore[i][j-i] );
		}

		if( nkozo )
		{
			ien = njob-1;
			ik = 0;
			for( i=0; i<ien; i++ )
			{
				jk = ik+1;
				for( j=i+1; j<njob; j++ ) 
				{
					if( kozoarivec[i] && kozoarivec[j] )
					{
						iscore_kozo[ik][jk-ik] = iscore[i][j-i];
					}
					if( kozoarivec[j] ) jk++;
				}
				if( kozoarivec[i] ) ik++;
			}
		}

//		fprintf( stderr, "Constructing a UPGMA tree ... " );
		fflush( stderr );
		if( topin )
		{
			fprintf( stderr, "--topin has been disabled\n" );
			exit( 1 );
//			fprintf( stderr, "Loading a topology ... " );
//			loadtop( njob, iscore, topol, len );
//			fprintf( stderr, "\ndone.\n\n" );
		}
		else if( subalignment ) // merge error no tame
		{
			fprintf( stderr, "Constructing a UPGMA tree ... " );
			fixed_supg_double_realloc_nobk_halfmtx_treeout_constrained( njob, iscore, topol, len, name, nlen, dep, nsubalignments, subtable, 1 );
		}
		else if( tbutree == 0 && compacttree ) // tbutree != 0 no toki (aln->mtx) ha, 6merdistance -> disttbfast.c; dp distance -> muzukashii
		{
			reporterr(       "Constructing a tree ... " );
			compacttree_memsaveselectable( njob, partmtx, mindistfrom, mindist, NULL, selfscore, seq, skiptable, topol, len, name, NULL, dep, treeout, compacttree, 1 );
			if( mindistfrom ) free( mindistfrom ); mindistfrom = NULL;
			if( mindist ) free( mindist );; mindist = NULL;
//			if( selfscore ) free( selfscore ); selfscore = NULL; // matomete free
			if( skiptable) FreeIntMtx( skiptable ); skiptable = NULL; // nikaime dake
			free( partmtx );
		}
		else if( treeout ) // merge error no tame
		{
			fprintf( stderr, "Constructing a UPGMA tree ... " );
			fixed_musclesupg_double_realloc_nobk_halfmtx_treeout( njob, iscore, topol, len, name, nlen, dep, 1 );
		}
		else
		{
			fprintf( stderr, "Constructing a UPGMA tree ... " );
			fixed_musclesupg_double_realloc_nobk_halfmtx( njob, iscore, topol, len, dep, 1, 1 );
		}
//		else 
//			ErrorExit( "Incorrect tree\n" );

		if( nkozo )
		{
//			for( i=0; i<nkozo-1; i++ )
//				for( j=i+1; j<nkozo; j++ )
//					fprintf( stderr, "iscore_kozo[%d][%d] =~ %f\n", i, j, iscore_kozo[i][j-i] );
			fixed_musclesupg_double_realloc_nobk_halfmtx( nkozo, iscore_kozo, topol_kozo, len_kozo, NULL, 1, 1 );
		}
		fprintf( stderr, "\ndone.\n\n" );
		fflush( stderr );
	}



	orderfp = fopen( "order", "w" );
	if( !orderfp )
	{
		fprintf( stderr, "Cannot open 'order'\n" );
		exit( 1 );
	}
	for( i=0; (j=topol[njob-2][0][i])!=-1; i++ )
	{
		fprintf( orderfp, "%d\n", j );
	}
	for( i=0; (j=topol[njob-2][1][i])!=-1; i++ )
	{
		fprintf( orderfp, "%d\n", j );
	}
	fclose( orderfp );

	if( treeout && noalign ) 
	{
		writeData_pointer( prep_g, njob, name, nlen, seq );
		fprintf( stderr, "\n" ); 
		SHOWVERSION;
		goto chudan; // 2016Jul31
	}

//	countnode( njob, topol, node0 );
	if( tbrweight )
	{
		weight = 3; 
#if 0
		utree = 0; counteff( njob, topol, len, eff ); utree = 1;
#else
		counteff_simple_double_nostatic( njob, topol, len, eff );
		for( i=njob-nadd; i<njob; i++ ) eff[i] /= (double)100;
#if 0
		fprintf( stderr, "######  WEIGHT = \n" );
		for( i=0; i<njob; i++ )
		{
			fprintf( stderr, "w[%d] = %f\n", i, eff[i] );
		}
		exit( 1 );
#endif
		if( nkozo )
		{
//			counteff_simple_double( nkozo, topol_kozo, len_kozo, eff_kozo ); // single weight nanode iranai
			for( i=0,j=0; i<njob; i++ )
			{
				if( kozoarivec[i] )
				{
//					eff_kozo_mapped[i] = eff_kozo[j]; //
					eff_kozo_mapped[i] = eff[i];      // single weight
					j++;
				}
				else
					eff_kozo_mapped[i] = 0.0;
//				fprintf( stderr, "eff_kozo_mapped[%d] = %f\n", i, eff_kozo_mapped[i] );
//				fprintf( stderr, "            eff[%d] = %f\n", i, eff[i] );
			}
		}


#endif
	}
	else
	{
		for( i=0; i<njob; i++ ) eff[i] = 1.0;
		if( nkozo ) 
		{
			for( i=0; i<njob; i++ ) 
			{
				if( kozoarivec[i] ) 
					eff_kozo_mapped[i] = 1.0;
				else
					eff_kozo_mapped[i] = 0.0;
			}
		}
	}

	if( iscore ) FreeFloatHalfMtx( iscore, njob ); iscore = NULL;
	FreeFloatMtx( len );

	alloclen = nlenmax*2+1; //chuui!
	bseq = AllocateCharMtx( njob, alloclen );


	if( nadd )
	{
		alignmentlength = strlen( seq[0] );
		for( i=0; i<njob-nadd; i++ )
		{
			if( alignmentlength != strlen( seq[i] ) )
			{
				fprintf( stderr, "#################################################################################\n" );
				fprintf( stderr, "# ERROR!                                                                        #\n" );
				fprintf( stderr, "# The original%4d sequences must be aligned                                    #\n", njob-nadd );
				fprintf( stderr, "#################################################################################\n" );
				exit( 1 );
			}
		}
		if( addprofile )
		{
			alignmentlength = strlen( seq[njob-nadd] );
			for( i=njob-nadd; i<njob; i++ )
			{
				if( alignmentlength != strlen( seq[i] ) )
				{
					fprintf( stderr, "###############################################################################\n" );
					fprintf( stderr, "# ERROR!                                                                      #\n" );
					fprintf( stderr, "# The%4d additional sequences must be aligned                                #\n", nadd );
					fprintf( stderr, "# Otherwise, try the '--add' option, instead of '--addprofile' option.        #\n" );
					fprintf( stderr, "###############################################################################\n" );
					exit( 1 );
				}
			}
			for( i=0; i<nadd; i++ ) addmem[i] = njob-nadd+i;
			addmem[nadd] = -1;
			foundthebranch = 0;
			for( i=0; i<njob-1; i++ )
			{
				if( samemember( topol[i][0], addmem ) ) // jissainiha nai
				{
					mergeoralign[i] = '1';
					foundthebranch = 1;
				}
				else if( samemember( topol[i][1], addmem ) ) // samemembern ni henkou kanou
				{
					mergeoralign[i] = '2';
					foundthebranch = 1;
				}
				else
				{
					mergeoralign[i] = 'n';
				}
			}
			if( !foundthebranch )
			{
				fprintf( stderr, "###############################################################################\n" );
				fprintf( stderr, "# ERROR!                                                                      #\n" );
				fprintf( stderr, "# There is no appropriate position to add the%4d sequences in the guide tree.#\n", nadd );
				fprintf( stderr, "# Check whether the%4d sequences form a monophyletic cluster.                #\n", nadd );
				fprintf( stderr, "# If not, try the '--add' option, instead of the '--addprofile' option.       #\n" );
				fprintf( stderr, "############################################################################### \n" );
				exit( 1 );
			}
			commongappick( nadd, seq+njob-nadd );
			for( i=njob-nadd; i<njob; i++ ) strcpy( bseq[i], seq[i] );
		}
		else
		{
			for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'n';
			for( j=njob-nadd; j<njob; j++ )
			{
				addmem[0] = j;
				addmem[1] = -1;
				for( i=0; i<njob-1; i++ )
				{
					if( samemembern( topol[i][0], addmem, 1 ) ) // arieru
					{
//						fprintf( stderr, "HIT!\n" );
						if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
						else mergeoralign[i] = '1';
					}
					else if( samemembern( topol[i][1], addmem, 1 ) )
					{
//						fprintf( stderr, "HIT!\n" );
						if( mergeoralign[i] != 'n' ) mergeoralign[i] = 'w';
						else mergeoralign[i] = '2';
					}
				}
			}
	
			for( i=0; i<nadd; i++ ) addmem[i] = njob-nadd+i;
			addmem[nadd] = -1;
			for( i=0; i<njob-1; i++ )
			{
				if( includemember( topol[i][0], addmem ) && includemember( topol[i][1], addmem ) )
				{
					mergeoralign[i] = 'w';
				}
				else if( includemember( topol[i][0], addmem ) )
				{
					mergeoralign[i] = '1';
				}
				else if( includemember( topol[i][1], addmem ) )
				{
					mergeoralign[i] = '2';
				}
			}
#if 0
			for( i=0; i<njob-1; i++ )
			{
				fprintf( stderr, "mem0 = " );
				for( j=0; topol[i][0][j]>-1; j++ )	fprintf( stderr, "%d ", topol[i][0][j] );
				fprintf( stderr, "\n" );
				fprintf( stderr, "mem1 = " );
				for( j=0; topol[i][1][j]>-1; j++ )	fprintf( stderr, "%d ", topol[i][1][j] );
				fprintf( stderr, "\n" );
				fprintf( stderr, "i=%d, mergeoralign[] = %c\n", i, mergeoralign[i] );
			}
#endif
			for( i=njob-nadd; i<njob; i++ ) gappick0( bseq[i], seq[i] );
		}

		commongappick( njob-nadd, seq );
		for( i=0; i<njob-nadd; i++ ) strcpy( bseq[i], seq[i] );
	}
//--------------- kokokara ----
	else if( subalignment )
	{
		for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'a';
		for( i=0; i<nsubalignments; i++ )
		{
			fprintf( stderr, "Checking subalignment %d:\n", i+1 );
			alignmentlength = strlen( seq[subtable[i][0]] );
//			for( j=0; subtable[i][j]!=-1; j++ )
//				fprintf( stderr, " %d. %-30.30s\n", subtable[i][j]+1, name[subtable[i][j]]+1 );
			for( j=0; subtable[i][j]!=-1; j++ )
			{
				if( subtable[i][j] >= njob )
				{
					fprintf( stderr, "No such sequence, %d.\n", subtable[i][j]+1 );
					exit( 1 );
				}
				if( alignmentlength != strlen( seq[subtable[i][j]] ) )
				{
					fprintf( stderr, "\n" );
					fprintf( stderr, "###############################################################################\n" );
					fprintf( stderr, "# ERROR!\n" );
					fprintf( stderr, "# Subalignment %d must be aligned.\n", i+1 );
					fprintf( stderr, "# Please check the alignment lengths of following sequences.\n" );
					fprintf( stderr, "#\n" );
					fprintf( stderr, "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][0]+1, name[subtable[i][0]]+1, alignmentlength );
					fprintf( stderr, "# %d. %-10.10s -> %d letters (including gaps)\n", subtable[i][j]+1, name[subtable[i][j]]+1, (int)strlen( seq[subtable[i][j]] ) );
					fprintf( stderr, "#\n" );
					fprintf( stderr, "# See http://mafft.cbrc.jp/alignment/software/merge.html for details.\n" );
					if( subalignmentoffset )
					{
						fprintf( stderr, "#\n" );
						fprintf( stderr, "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
						fprintf( stderr, "# In this case, the rule of numbering is:\n" );
						fprintf( stderr, "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
						fprintf( stderr, "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
					}
					fprintf( stderr, "###############################################################################\n" );
					fprintf( stderr, "\n" );
					exit( 1 );
				}
				insubtable[subtable[i][j]] = 1;
			}
			for( j=0; j<njob-1; j++ )
			{
				if( includemember( topol[j][0], subtable[i] ) && includemember( topol[j][1], subtable[i] ) )
				{
					mergeoralign[j] = 'n';
				}
			}
			foundthebranch = 0;
			for( j=0; j<njob-1; j++ )
			{
				if( samemember( topol[j][0], subtable[i] ) || samemember( topol[j][1], subtable[i] ) )
				{
					foundthebranch = 1;
					fprintf( stderr, " -> OK\n" );
					break;
				}
			}
			if( !foundthebranch )
			{
				system( "cp infile.tree GuideTree" ); // tekitou
				fprintf( stderr, "\n" );
				fprintf( stderr, "###############################################################################\n" );
				fprintf( stderr, "# ERROR!\n" );
				fprintf( stderr, "# Subalignment %d does not form a monophyletic cluster\n", i+1 );
				fprintf( stderr, "# in the guide tree ('GuideTree' in this directory) internally computed.\n" );
				fprintf( stderr, "# If you really want to use this subalignment, pelase give a tree with --treein \n" );
				fprintf( stderr, "# http://mafft.cbrc.jp/alignment/software/treein.html\n" );
				fprintf( stderr, "# http://mafft.cbrc.jp/alignment/software/merge.html\n" );
				if( subalignmentoffset )
				{
					fprintf( stderr, "#\n" );
					fprintf( stderr, "# You specified seed alignment(s) consisting of %d sequences.\n", subalignmentoffset );
					fprintf( stderr, "# In this case, the rule of numbering is:\n" );
					fprintf( stderr, "#   The aligned seed sequences are numbered as 1 .. %d\n", subalignmentoffset );
					fprintf( stderr, "#   The input sequences to be aligned are numbered as %d .. %d\n", subalignmentoffset+1, subalignmentoffset+njob );
				}
				fprintf( stderr, "############################################################################### \n" );
				fprintf( stderr, "\n" );
				exit( 1 );
			}
//			commongappick( seq[subtable[i]], subalignment[i] ); // irukamo
		}
#if 0
		for( i=0; i<njob-1; i++ )
		{
			fprintf( stderr, "STEP %d\n", i+1 );
			fprintf( stderr, "group1 = " );
			for( j=0; topol[i][0][j] != -1; j++ )
				fprintf( stderr, "%d ", topol[i][0][j]+1 );
			fprintf( stderr, "\n" );
			fprintf( stderr, "group2 = " );
			for( j=0; topol[i][1][j] != -1; j++ )
				fprintf( stderr, "%d ", topol[i][1][j]+1 );
			fprintf( stderr, "\n" );
			fprintf( stderr, "%d -> %c\n\n", i, mergeoralign[i] );
		}
#endif

		for( i=0; i<njob; i++ ) 
		{
			if( insubtable[i] ) strcpy( bseq[i], seq[i] );
			else gappick0( bseq[i], seq[i] );
		}

		for( i=0; i<nsubalignments; i++ ) 
		{
			for( j=0; subtable[i][j]!=-1; j++ ) subalnpt[i][j] = bseq[subtable[i][j]];
			if( !preservegaps[i] ) commongappick( j, subalnpt[i] );
		}

		FreeIntMtx( subtable );
		free( insubtable );
		for( i=0; i<nsubalignments; i++ ) free( subalnpt[i] );
		free( subalnpt );
		free( preservegaps );
	}
//--------------- kokomade ----
	else
	{
		for( i=0; i<njob; i++ ) gappick0( bseq[i], seq[i] );
		for( i=0; i<njob-1; i++ ) mergeoralign[i] = 'a';
	}

	if( rnakozo && rnaprediction == 'm' )
	{
		singlerna = (RNApair ***)calloc( njob, sizeof( RNApair ** ) );
		prep = fopen( "hat4", "r" );
		if( prep == NULL ) ErrorExit( "Make hat4 using mccaskill." );
		fprintf( stderr, "Loading 'hat4' ... " );
		for( i=0; i<njob; i++ )
		{
			nogaplen = strlen( bseq[i] );
			singlerna[i] = (RNApair **)calloc( nogaplen+1, sizeof( RNApair * ) );
			for( j=0; j<nogaplen; j++ )
			{
				singlerna[i][j] = (RNApair *)calloc( 1, sizeof( RNApair ) );
				singlerna[i][j][0].bestpos = -1;
				singlerna[i][j][0].bestscore = -1.0;
			}
			singlerna[i][nogaplen] =  NULL;
//			fprintf( stderr, "### reading bpp %d ...\n", i );
			readmccaskill( prep, singlerna[i], nogaplen );
		}
		fclose( prep );
		fprintf( stderr, "\ndone.\n" );
	}
	else
		singlerna = NULL;


	fprintf( stderr, "Progressive alignment ... \n" );

#ifdef enablemultithread
	if( nthread > 0 && nadd == 0 )
	{
		treebasethread_arg_t *targ;	
		int jobpos;
		pthread_t *handle;
		pthread_mutex_t mutex;
		pthread_cond_t treecond;
		int *fftlog;
		int nrun;
		int nthread_yoyu;

		nthread_yoyu = nthread * 1;
		nrun = 0;
		jobpos = 0;
		targ = calloc( nthread_yoyu, sizeof( treebasethread_arg_t ) );
		fftlog = AllocateIntVec( njob );
		handle = calloc( nthread_yoyu, sizeof( pthread_t ) );
		pthread_mutex_init( &mutex, NULL );
		pthread_cond_init( &treecond, NULL );

		for( i=0; i<njob; i++ ) dep[i].done = 0;
		for( i=0; i<njob; i++ ) fftlog[i] = 1;

		if( constraint )
		{
			if( specifictarget )
				calcimportance_target( njob, ntarget, eff, bseq, localhomtable, targetmap, targetmapr );
//				dontcalcimportance_target( njob, eff, bseq, localhomtable, ntarget ); // CHUUIII!!!!!
			else
				calcimportance_half( njob, eff, bseq, localhomtable );
		}
//			dontcalcimportance( njob, eff, bseq, localhomtable ); // CHUUUUIIII!!!

		for( i=0; i<nthread_yoyu; i++ )
		{
			targ[i].thread_no = i;
			targ[i].nrunpt = &nrun;
			targ[i].njob = njob;
			targ[i].nlen = nlen;
			targ[i].jobpospt = &jobpos;
			targ[i].topol = topol;
			targ[i].dep = dep;
			targ[i].aseq = bseq;
			targ[i].effarr = eff;
			targ[i].alloclenpt = &alloclen;
			targ[i].localhomtable = localhomtable;
			targ[i].singlerna = singlerna;
			targ[i].effarr_kozo = eff_kozo_mapped;
			targ[i].fftlog = fftlog;
			targ[i].mergeoralign = mergeoralign;
			targ[i].targetmap = targetmap;
			targ[i].mutex = &mutex;
			targ[i].treecond = &treecond;

			pthread_create( handle+i, NULL, treebasethread, (void *)(targ+i) );
		}

		for( i=0; i<nthread_yoyu; i++ )
		{
			pthread_join( handle[i], NULL );
		}
		pthread_mutex_destroy( &mutex );
		pthread_cond_destroy( &treecond );
		free( handle );
		free( targ );
		free( fftlog );
//		free( topol[njob-1][0] );
//		free( topol[njob-1][1] );
		free( topol[njob-1] );
		free( topol );
	}
	else
#endif

		treebase( nlen, bseq, nadd, mergeoralign, mseq1, mseq2, topol, dep, eff, &alloclen, localhomtable, singlerna, eff_kozo_mapped, targetmap, targetmapr, ntarget );
	fprintf( stderr, "\ndone.\n" );

	
	if( keeplength )
	{

		dlf = fopen( "_deletelist", "w" );
		deletelist = (int **)calloc( nadd+1, sizeof( int * ) );
		for( i=0; i<nadd; i++ )
		{
			deletelist[i] = calloc( 1, sizeof( int ) );
			deletelist[i][0] = -1;
		}
		deletelist[nadd] = NULL;
		ndeleted = deletenewinsertions_whole( njob-nadd, nadd, bseq, bseq+njob-nadd, deletelist );

		for( i=0; i<nadd; i++ )
		{
			if( deletelist[i] )
				for( j=0; deletelist[i][j]!=-1; j++ )
					fprintf( dlf, "%d %d\n", njob-nadd+i, deletelist[i][j] ); // 0origin
		}
		fclose( dlf );

		restoreoriginalgaps( njob, bseq, originalgaps );

		if( mapout )
		{
			dlf = fopen( "_deletemap", "w" );
			reconstructdeletemap( nadd, addbk, deletelist, bseq+njob-nadd, dlf, name+njob-nadd );
			FreeCharMtx( addbk );
			addbk = NULL;
			fclose( dlf );
		}

		FreeIntMtx( deletelist );
		deletelist = NULL;
	}



	if( scoreout )
	{
		unweightedspscore = plainscore( njob, bseq );
		fprintf( stderr, "\nSCORE %s = %.0f, ", "(treebase)", unweightedspscore );
		fprintf( stderr, "SCORE / residue = %f", unweightedspscore / ( njob * strlen( bseq[0] ) ) );
		fprintf( stderr, "\n\n" );
	}

#if 0
	if( constraint )
	{
		LocalHom *tmppt1, *tmppt2;
		for( i=0; i<njob; i++ )
		{
			for( j=0; j<njob; j++ )
			{
				tmppt1 = localhomtable[i]+j;
				while( tmppt2 = tmppt1->next )
				{
					free( (void *)tmppt1 );
					tmppt1 = tmppt2;
				}
				free( (void *)tmppt1 );
			}
			free( (void *)(localhomtable[i]+j) );
		}
		free( (void *)localhomtable );
	}
#endif

	fprintf( trap_g, "done.\n" );
//	fclose( trap_g );
	free( mergeoralign );
	freeconstants();



	if( rnakozo && rnaprediction == 'm' ) 
	{
		if( singlerna ) // nen no tame
		{
			for( i=0; i<njob; i++ ) 
			{
				for( j=0; singlerna[i][j]!=NULL; j++ )
				{
					if( singlerna[i][j] ) free( singlerna[i][j] );
				}
				if( singlerna[i] ) free( singlerna[i] );
			}
			free( singlerna );
			singlerna = NULL;
		}
	}

	writeData_pointer( prep_g, njob, name, nlen, bseq );
#if 0
	writeData( stdout, njob, name, nlen, bseq );
	writePre( njob, name, nlen, bseq, !contin );
	writeData_pointer( prep_g, njob, name, nlen, aseq );
#endif
#if IODEBUG
	fprintf( stderr, "OSHIMAI\n" );
#endif

	if( constraint ) 
	{
		if( specifictarget )
			FreeLocalHomTable_part( localhomtable, ntarget, njob );
		else
			FreeLocalHomTable_half( localhomtable, njob );
	}
	free( targetmap );
	free( targetmapr );


	if( spscoreout ) reporterr( "Unweighted sum-of-pairs score = %10.5f\n", sumofpairsscore( njob, bseq ) );
	SHOWVERSION;
	if( ndeleted > 0 )
	{
		reporterr( "\nTo keep the alignment length, %d letters were DELETED.\n", ndeleted );
		if( mapout )
			reporterr( "The deleted letters are shown in the (filename).map file.\n" );
		else
			reporterr( "To know the positions of deleted letters, rerun the same command with the --mapout option.\n" );
	}


	free( kozoarivec );
	FreeCharMtx( seq );
	FreeCharMtx( bseq );
	free( mseq1 );
	free( mseq2 );

	FreeCharMtx( name );
	free( nlen );
	free( selfscore );

//	for( i=0; i<njob; i++ ) 
//	{
//		free( topol[i][0] );
//		free( topol[i][1] );
//		free( topol[i] );
//	}
//	free( topol );
//	free( len );
//	free( iscore );
	free( eff );
	free( dep );
	closeFiles();
	if( nadd ) free( addmem );

	return( 0 );

chudan:
	if( seq ) FreeCharMtx( seq ); seq = NULL;
	if( mseq1 ) free( mseq1 ); mseq1 = NULL;
	if( mseq2 ) free( mseq2 ); mseq2 = NULL;

	if( name ) FreeCharMtx( name ); name = NULL;
	if( nlen ) free( nlen ); nlen = NULL;
	if( selfscore ) free( selfscore ); selfscore = NULL;
	if( mergeoralign ) free( mergeoralign ); mergeoralign = NULL;

	if( localhomtable )
	{
		reporterr( "freeing localhomtable\n" );
		if( specifictarget )
			FreeLocalHomTable_part( localhomtable, ntarget, njob );
		else
			FreeLocalHomTable_half( localhomtable, njob );
	} 
	localhomtable = NULL;
	if( targetmap ) free( targetmap ); targetmap = NULL;
	if( targetmapr ) free( targetmapr ); targetmapr = NULL;

	if( kozoarivec ) free( kozoarivec ); kozoarivec = NULL;


	if( topol ) FreeIntCub( topol ); topol = NULL;
	if( len ) FreeFloatMtx( len ); len = NULL;
	if( iscore ) FreeFloatHalfMtx( iscore, njob ); iscore = NULL;
	if( eff ) free( eff ); eff = NULL;
	if( dep ) free( dep ); dep = NULL;

	freeconstants();
	closeFiles();
	FreeCommonIP();
	return( 0 );

}
