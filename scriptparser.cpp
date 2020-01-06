#include "scriptparser.h"
#include "vstring.h"

#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdlib>
#include <cassert>
#include <list>

using namespace std;


void ScriptParser::getStrPairs(const char *srcfn, vector<const char*> &strlst1,
                                     vector<const char*> &strlst2)
{
    assert(srcfn);
    printf("Read string pairs ........... ");
	FILE *fp = fopen(srcfn, "r");
	if(fp == NULL)
	{
		printf("File '%s' cannot be open!\n", srcfn);
		return ;
	}

	char result[STRLEN];
	char former[STRLEN];
	char latter[STRLEN];

	char *first = NULL, *second = NULL;
	unsigned int i = 0, len1 = 0, len2 = 0;
	int nargs = 0;

	while(!feof(fp))
	{
		nargs = fscanf(fp, "%s\n", result);

		if(nargs == 0)
		continue;

		if((!strcmp(result, "NULL"))||(!strcmp(result, "")))
		{
		      continue;
		}
		else
		{
		    VString::split_twin(former, latter, result, '^');
		    len1   = strlen(former) + 1;
		    len2   = strlen(latter) + 1;
		    first  = new char[len1];
		    second = new char[len2];

		    strcpy(first,  former);
		    strcpy(second, latter);

		    strlst1.push_back(first);
		    strlst2.push_back(second);
		    i++;
		}
		strcpy(result, "");
	}
	printf("%d\n", i);
	fclose(fp);
	return ;
}

vector<const char*> ScriptParser::getNamelst(const char *fname)
{
    assert(fname);
    vector<const char*> kflst;

	FILE *fp = fopen(fname, "r");
	if(fp == NULL)
	{
		printf("File '%s' cannot be open!\n",fname);
		return kflst;
	}

	char result[ScriptParser::STRLEN];
	char *crntKeyFrm = NULL;
	int i = 0, len = 0, nargs = 0;

	while(!feof(fp))
	{
		nargs = fscanf(fp, "%s\n", result);

		if(nargs == 0)
		continue;

		len = strlen(result) + 1;
		if((!strcmp(result, "NULL"))||(!strcmp(result, "")))
		{
		      continue;
		}
		else
		{
             crntKeyFrm = new char[len];
             strcpy(crntKeyFrm, result);
		     kflst.push_back(crntKeyFrm);
		     i++;
		}
		strcpy(result, "");
	}
	fclose(fp);
	return kflst;
}


void ScriptParser::test()
{
	vector<const char *> kflst1;
	vector<const char *> kflst2;

	const char *fn    = "e:/wanlei/vcdlab/itmtabs.txt";

	map<unsigned int, map<string, const char*>* > itmMaps;
	itmMaps = ScriptParser::getItmMaps(fn);

	map<unsigned int, map<string, const char*>* >::iterator it;
	map<string, const char*>::iterator iter;
	map<string, const char*> *crntMap;
	string keystr;
	for(it = itmMaps.begin(); it != itmMaps.end(); it++)
	{
	    cout<<it->first<<endl;
	    crntMap = it->second;
	    for(iter = crntMap->begin(); iter != crntMap->end(); iter++)
	    {
	        cout<<iter->first<<"\t"<<iter->second<<endl;
	        keystr =  iter->first;
	        delete [] iter->second;
	        keystr.clear();
	    }
	    cout<<endl;
	    crntMap->clear();
	    delete crntMap;
	}
	itmMaps.clear();
}
