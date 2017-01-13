#include "mltaln.h"
int main( int argc, char **argv ) //may be this is the main method called for pairwise alignment. It calls 'pairlocalalign' method in pairlocalalign.c file
{
	int res = pairlocalalign( 0, 0, NULL, NULL, NULL, NULL, argc, argv );
	if( res == GUI_CANCEL ) res = 0; // treeout de goto chudan wo riyousuru
	return res;
}
