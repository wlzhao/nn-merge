#ifndef SCRIPTPARSER_H
#define SCRIPTPARSER_H

#include "vstring.h"

#include <cstring>
#include <string>
#include <vector>
#include <cassert>
#include <list>
#include <map>

using namespace std;


class ScriptParser
{

public:
    static const int STRLEN = 512;
	ScriptParser(){}
	static vector<const char*> getNamelst(const char *fname);
	static void getStrPairs(const char* srcfn, vector<const char*> &strlst1,
                            vector<const char*> &strlst2);

static map<string, const char *> getConf(const char *fn)
{
    assert(fn);
    FILE *fp = fopen(fn, "r");
	if(fp == NULL)
	{
		printf("Configure '%s' cannot be open!\n", fn);
		exit(0);
	}

	map<string, const char *> paras;

    string fst_str = "", sec_str = "";
    char *pairf = NULL, *pairs = NULL, result[STRLEN];
    char *ptStr = NULL;

	while(!feof(fp))
	{
		ptStr = fgets(result, STRLEN, fp);

		if(ptStr == NULL)
		continue;

		if(!strcmp(ptStr, ""))
		continue;

		if(result[0] == '#')
		continue;

		VString::trimEnd(result, '\n');
		VString::trimEnd(result, '\r');
		VString::parsePair(result, fst_str, sec_str, '=');

		pairf = new char[fst_str.length()+1];
		strcpy(pairf, fst_str.c_str());

		pairs = new char[sec_str.length()+1];
		strcpy(pairs, sec_str.c_str());

		paras.insert(pair<const char *, const char*>(pairf, pairs));

		fst_str.erase(fst_str.begin(), fst_str.end());
		sec_str.erase(sec_str.begin(), sec_str.end());
		strcpy(result, "");
	}
	fclose(fp);

	return paras;
}

static map<unsigned int, map<string,const char*>* > getItmMaps(const char *fn)
{
  assert(fn);

    unsigned int index = 0;
    FILE *fp = fopen(fn, "r");
	if(fp == NULL)
	{
		printf("Configure '%s' cannot be open!\n", fn);
		exit(0);
	}

	char result[STRLEN];

	map<unsigned int, map<string, const char*>* > data;
	map<string, const char*> *crnt_para_map;
    string fst_str = "", sec_str = "";
    char *pairf = NULL, *pairs = NULL;
    char *ptStr = NULL;

	while(!feof(fp))
	{
		ptStr = fgets(result, STRLEN, fp);

		if(ptStr == NULL)
		continue;

		if(!strcmp(ptStr, ""))
		continue;

		VString::trimEnd(result, '\n');
        VString::trimEnd(result, '\r');

		if(!strcmp(result, "<item>"))
		{
		    crnt_para_map = new map<string, const char*>;
		    ptStr = fgets(result, STRLEN, fp);

		    if(!strcmp(result, ""))
		    continue;

		    if(result[0] == '#')
		    continue;

		    VString::trimEnd(result, '\n');
		    VString::trimEnd(result, '\r');
		    while(strcmp(result, "</item>"))
		    {
		        VString::parsePair(result, fst_str, sec_str, '=');
		        pairf = new char[fst_str.length()+1];
		        strcpy(pairf, fst_str.c_str());

		        pairs = new char[sec_str.length()+1];
		        strcpy(pairs, sec_str.c_str());
		        crnt_para_map->insert(pair<const char*,const char*>(pairf,pairs));
		        fst_str.erase(fst_str.begin(),fst_str.end());
		        sec_str.erase(sec_str.begin(),sec_str.end());

		        ptStr = fgets(result, STRLEN, fp);

		        if(!strcmp(ptStr, ""))
		        continue;

		        VString::trimEnd(result, '\n');
		        VString::trimEnd(result, '\r');
		    }
		    data.insert(pair<int,map<string,const char*>* >(index, crnt_para_map));
		    index++;
		}
		strcpy(result, "");
	}
	fclose(fp);
	return data;
}

static void clearParaMap(map<string, const char *> &paras)
{
    string key;
    const char *val;
    map<string, const char *>::iterator mit;

    for(mit = paras.begin(); mit != paras.end(); mit++)
    {
        key = mit->first;
        key.clear();
        val = mit->second;
        delete [] val;
    }
    paras.clear();
}

static void clearItmMaps(map<unsigned int, map<string, const char*>* > &itmMaps)
{
    unsigned int idx = 0;
    map<string, const char*>* crntMap;
    for(idx  = 0; idx < itmMaps.size(); idx++)
    {
        crntMap = itmMaps[idx];
        ScriptParser::clearParaMap(*crntMap);
        delete crntMap;
    }
    itmMaps.clear();
}

static bool verifyParaMap(const char *options[], const unsigned int num,
                               map<string, const char *> &paras)
{
    unsigned int i;
    for(i = 0; i < num; i++)
    {
        if(paras.find(options[i]) == paras.end())
        {
            cout<<"Option '"<<options[i]<<"' has not been set!\n";
        }
    }
    return true;
}
	virtual ~ScriptParser(){}
	static void test();


};

#endif
