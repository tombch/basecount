#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <vector>
#include <string>


namespace py = pybind11;


std::vector<std::vector<unsigned int>> compute(unsigned int refLen, std::vector<std::string> reads, std::vector<unsigned int> starts, std::vector<std::vector<std::pair<unsigned int, unsigned int>>> ctuples)
{
    std::vector<std::vector<unsigned int>> baseCounts(refLen, std::vector<unsigned int>(6));

    for (unsigned int i = 0; i < reads.size(); i++)
    {
        unsigned int seqPos = 0;
        unsigned int refPos = starts[i];
        auto &read = reads[i];
        auto &tups = ctuples[i];

        for (auto &tup : tups)
        {
            auto &operation = tup.first;
            auto &opLen = tup.second;
            if ((operation == 0) || (operation == 7) || (operation == 8))
            {
                for (unsigned int j = 0; j < opLen; j++)
                {
                    switch(read[seqPos + j])
                    {
                        case 'A' : baseCounts.at(refPos + j).at(0) += 1; break;
                        case 'C' : baseCounts.at(refPos + j).at(1) += 1; break;
                        case 'G' : baseCounts.at(refPos + j).at(2) += 1; break;
                        case 'T' : baseCounts.at(refPos + j).at(3) += 1; break;
                        case 'N' : baseCounts.at(refPos + j).at(4) += 1; break;
                    }
                }

                refPos += opLen;
                seqPos += opLen;
            }
            else if (operation == 1)
                seqPos += opLen;

            else if ((operation == 2) || (operation == 3))
            {
                for (unsigned int j = 0; j < opLen; j++)
                    baseCounts.at(refPos + j).at(5) += 1;

                refPos += opLen;
            }
        }
    }
    return baseCounts;
}


PYBIND11_MODULE(basecomp, m)
{
    m.doc() = "Count bases";
    m.def("compute", &compute, "Count bases");
}
