package com.hartwig.hmftools.geneutils.ensembl;

import static com.hartwig.hmftools.common.utils.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.geneutils.common.CommonUtils.GU_LOGGER;

import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Maps;

public class  HgncGenes
{
    private final Map<String,HgncGene> mDataByHgncId; // keyed by HGNC Id
    private final Map<String,HgncGene> mDataByGeneId;
    private final Map<String,HgncGene> mDataBySymbol;

    public HgncGenes(final String filename)
    {
        mDataByHgncId = Maps.newHashMap();
        mDataByGeneId = Maps.newHashMap();
        mDataBySymbol = Maps.newHashMap();

        loadGeneData(filename);
    }

    public boolean hasData() { return !mDataByHgncId.isEmpty(); }

    public HgncGene getByHgncId(final String hgncId)
    {
        return mDataByHgncId.get(hgncId);
    }

    public HgncGene getByGeneId(final String geneId)
    {
        return mDataByGeneId.get(geneId);
    }

    public HgncGene getBySymbol(final String symbol)
    {
        return mDataBySymbol.get(symbol);
    }

    private void loadGeneData(final String filename)
    {
        if(filename == null)
            return;

        try
        {
            List<String> lines = Files.readAllLines(Paths.get(filename));

            final Map<String, Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), "\t");
            lines.remove(0);

            int hgncIdIndex = fieldsIndexMap.get("hgnc_id");
            int ensemblGeneIdIndex = fieldsIndexMap.get("ensembl_gene_id");
            int symbolIndex = fieldsIndexMap.get("symbol");

            for(String line : lines)
            {
                String[] values = line.split("\t", -1);

                if(values.length < 20)
                    continue;

                String hgncId = values[hgncIdIndex];
                String ensemblGeneId = values[ensemblGeneIdIndex];
                String symbol = values[symbolIndex];

                if(ensemblGeneId.isEmpty())
                    continue;

                HgncGene geneData = new HgncGene(hgncId, ensemblGeneId, symbol);

                mDataByHgncId.put(hgncId, geneData);
                mDataByGeneId.put(ensemblGeneId, geneData);
                mDataBySymbol.put(symbol, geneData);
            }

            GU_LOGGER.info("loaded {} HGNC genes from file({})", mDataByHgncId.size(), filename);
        }
        catch(Exception e)
        {
            GU_LOGGER.error("failed to read HGNC gene file({}): {}", filename, e.toString());
        }
    }

}
