//
//  main.cpp
//  설계
//
//  Created by 정용헌 on 2020/05/28.
//  Copyright © 2020 정용헌. All rights reserved.
//
// 설계 프로젝트 - 2015112113 Algorithm

#include <iostream>
#include <string>
#include <chrono>
#include <math.h>
#include <fstream>
#include <vector>

using namespace std;

#define N 500000 // up to 3000000000
#define M 50000 // up to 200milion ~ 2000 milion
#define L 30 // 30 ~ 60
#define D 3 // Snip

ofstream output("DNAinfo.txt"); // txt파일로 output
chrono::system_clock::time_point StartTime; // 실행시간 계산을 위한 chrono
chrono::system_clock::time_point EndTime;

const int divN = (int)(L / (D + 1)); // shortread를 나누는 길이 - 36 : 3일 경우 36/ 4 인 9


void getTime() // 실행시간을 구하기 위한 함수 설정
{
   chrono::milliseconds mill = chrono::duration_cast<chrono::milliseconds>(EndTime - StartTime);
   cout << "하는데 걸린 시간 : " << mill.count() << " ms" << endl;
}

int missNum(string& A, string& B) // mismatch  확인
{
   int miss = 0;
   for (int i = 0; i < N; i++)
   {
      if (A[i] != B[i])
         miss++;
   }
   return miss;
}


class DNAProj
{
private:

   string mySequence = " "; // 염기서열 Sequence 만들 string
   string referSequence = " ";// 염기서열 Sequecne 만들 string
   vector<string> shortRead;
   vector<vector<int>> indexTable;
   int ranBound = sqrt(N); // 난수 범위보다 늘어나는 N때문에 제곱하기 위해 계산

public:
    
   DNAProj()
    {
      StartTime = chrono::system_clock::now();
      shortRead.resize(M);
      makeGenome();
      makeRefer();
      EndTime = chrono::system_clock::now();
      cout << "DNA sequence 생성";
      getTime();

      StartTime = chrono::system_clock::now();
      makeShort();
      EndTime = chrono::system_clock::now();
      cout << "\nshort read 생성";
      getTime();
   }

   void getString(string& A, string& B, vector<string>& C) // 각자 dna 생성을 위한 스트링 생성
   {
      A = mySequence;
      B = referSequence;
      C = shortRead;
   }

   char AdaptC(int& input)
   {
      switch (input)
      {
      case 0:
         return 'A';
      case 1:
         return 'C';
      case 2:
         return 'G';
      case 3:
         return 'T';
      }
       return input;
   }

   int AdaptI(char& input)
   {
      switch (input)
      {
      case 'A':
         return 0;
      case 'C':
         return 1;
      case 'G':
         return 2;
      case 'T':
         return 3;
      }
       return input;
   }

   void makeGenome()
   {
      int ran;
      char temp;
      while (mySequence.size() <= N + 1) //사이즈 N만큼 염기서열 생성
      {
         ran = rand() % 4;
         temp = AdaptC(ran);
         mySequence += temp;
         referSequence += temp; // 생성되는 난수에 따라 뒤에 텍스트를 붙여줌
      }
   }
    
   void makeRefer() // reference DNA 생성
   {

      int ran1, ran2;
      for (int i = 0; i < (N * 0.05); i++)
      {
         ran1 = rand() % 4;
         ran2 = ((rand() % ranBound) + 1) * ((rand() % ranBound) + 1);
         referSequence[ran2] = AdaptC(ran1);
      }
      cout << "my DNA와 reference DNA 차이 : " << missNum(mySequence, referSequence)
         << " 개" << endl << endl;

      cout << "my DNA와 reference DNA 일치율 : " << (1 - (double)missNum(mySequence, referSequence) / N) * 100
         << "%" << endl << endl;
   }

   void makeShort() // shortread 생성
   {
      int ran;
      for (int i = 0; i < M; i++) {
         ran = ((rand() % ranBound) + 1) * ((rand() % ranBound) + 1); // 임의의 위치 정하기위해 난수 설정
   
         if (ran > N - L)
            ran -= L; // 배열 끝이 시퀀스를 넘어가면 안되므로 L을 빼줌
         
         for (int j = 0; j < L; j++)
            shortRead[i] += mySequence[ran++];
      }
   }
   void makeTable()
   {
      StartTime = chrono::system_clock::now();
      indexTable.resize((int)pow(4, divN));
      const int d = 4;
      const int m = divN; // shortread의 길이 : (L / (D + 1))
      int t = 0;
      int Dl = (int)(pow(d, m - 1));
      for (int i = 0; i < m; i++)
      {
         t = (d * t + AdaptI(referSequence[i + 1])); // 해시값 계산
      }
      indexTable[t].push_back(1); // push_back 함수 index 값 증가
      for (int i = 2; i < N - m + 2; i++) // 처음부터 끝까지 탐색
      {
         t = ((d * (t - AdaptI(referSequence[i - 1]) * Dl)) + (AdaptI(referSequence[i + m - 1])));
         indexTable[t].push_back(i);
      }   // 그다음 해쉬값 계산
      EndTime = chrono::system_clock::now();
      cout << "\nindex Table 생성";
      getTime();
   }

   void makeResult()
   {
      StartTime = chrono::system_clock::now();
      int index = 0;
      for (int i = 0; i < M; i++)
      {
         index = findshortIndex(i);
         if (index != -1)
         {
            for (int j = 0; j < L; j++)
            {
               referSequence[index++] = shortRead[i][j];
            }
         }

      }
      EndTime = chrono::system_clock::now();
      cout << "\n결과 생성";
      getTime();

      cout << "\n정확성 : " << (1 - (double)missNum(mySequence, referSequence) / N) * 100
         << "%" << endl << endl;
   }

   int findshortIndex(int input)
   {
      int t = 0;
      for (int i = 0; i < D + 1; i++) // D+1 = 4
      {
         t = 0;
         for (int j = 0; j < divN; j++) // divN = 9
         {
            t = (4 * t + AdaptI(shortRead[input][divN * i + j])); // short 4개씩 나눠서 라빈 카프
         }
         for (int k = 0; k < indexTable[t].size(); k++)
         {
            if (compareString(indexTable[t][k], input, i))
               return (indexTable[t][k] - divN * i);
         }
      }
      return -1;
   }

   bool compareString(int input, int index, int order)
   {
      int misMatch = 0;
      for (int i = 0; i < L; i++)
      {
         if (referSequence[input++ - divN * order] != shortRead[index][i])
            misMatch++;
         if (misMatch >= 4)
            return false;
      }
      return true;
   }

   void printInfo() // txt 파일에 생성된 랜덤 염기와 해시값, 그에 해당하는 주소를 출력한다.
   {
      for (int i = 1; i < N + 1; i++)
         output << referSequence[i];
      output << endl << endl;

      output << "my : ";
      for (int i = 1; i < N + 1; i++)
         output << mySequence[i];
      output << endl << endl;

      output << "shortread : ";
      for (int i = 0; i < M; i++) {
         output << i << " : ";
         for (int j = 0; j < L; j++)
         {
            output << shortRead[i][j];
         }
         output << endl;
      }
      output << endl << endl;

      output << "indexTable : ";
      for (int i = 0; i < (int)pow(4, divN); i++) {
         output << i << " : ";
         for (int j = 0; j < indexTable[i].size(); j++)
         {
            output << indexTable[i][j] << " ";
         }
         output << endl;
      }
      output << endl << endl;
   }
};

class bruteforce // 알고리즘 비교를 위한 trivial
{
private:
   string my;
   string refer;
   vector<string> shortR;

public:
   bruteforce(DNAProj G)
   {
      G.getString(my, refer, shortR);
   }

   void makeResult()
   {
      cout << "my와 refer 차이 : " << missNum(my, refer)
         << " 개" << endl << endl;

      StartTime = chrono::system_clock::now();
      int index = 0;
      for (int i = 0; i < M; i++)
      {
         index = findIndex(i);

         if (index != -1)
         {
            for (int j = 0; j < L; j++)
            {
               refer[index++] = shortR[i][j];
            }
         }
      }
      EndTime = chrono::system_clock::now();
      cout << "결과 생성";
      getTime();

      cout << "trivial 정확성 : " << (1 - (double)missNum(my, refer) / N) * 100
         << "%" << endl << endl;
   }

   int findIndex(int shortIndex)
   {
      int temp, miss;
      bool match = false;
      int index = 0;
      for (int i = 1; i < N - L + 1; i++)
      {
         miss = 0;
         index = i;
         temp = i;
         for (int j = 0; j < L; j++)
         {
            if (shortR[shortIndex][j] != refer[index++]) miss++;
            if (miss >= 4) break;
            if (j == L - 1) match = true;
         }
         if (match) return index - L;
      }
      return -1;
   }
};

int main()
{
   srand((unsigned int)time(NULL)); // 랜덤함수 중복방지 위함
   cout << endl << endl << " ====== 2015112113 정용헌 설계 프로젝트 =======" <<endl <<endl ;
   cout << " Rabin-Karp Algorithm" << endl << endl;
   DNAProj DNAProj;
   bruteforce brute(DNAProj);

   DNAProj.makeTable();
   DNAProj.printInfo();
   output.close();
   DNAProj.makeResult();
   cout << "=============================================" << endl << endl;
   cout << " trivial algorithm" << endl << endl;
   brute.makeResult();
   output.close();
   return 0;
}


