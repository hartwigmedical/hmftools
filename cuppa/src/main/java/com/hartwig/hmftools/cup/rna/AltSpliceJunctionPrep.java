package com.hartwig.hmftools.cup.rna;

import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_END;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_POS_START;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.cup.common.CupConstants.CUP_LOGGER;
import static com.hartwig.hmftools.cup.prep.DataSource.RNA;
import static com.hartwig.hmftools.cup.prep.PrepConfig.REF_ALT_SJ_SITES;

import java.io.BufferedReader;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.rna.NovelSpliceJunction;
import com.hartwig.hmftools.common.rna.NovelSpliceJunctionFile;
import com.hartwig.hmftools.cup.prep.CategoryType;
import com.hartwig.hmftools.cup.prep.CategoryPrep;
import com.hartwig.hmftools.cup.prep.DataItem;
import com.hartwig.hmftools.cup.prep.ItemType;
import com.hartwig.hmftools.cup.prep.PrepConfig;

public class AltSpliceJunctionPrep implements CategoryPrep
{
    private final PrepConfig mConfig;

    private final Map<String,Integer> mRefAsjIndexMap; // map from Alt-SJ into matrix rows, and which sites to use

    public AltSpliceJunctionPrep(final PrepConfig config)
    {
        mConfig = config;

        mRefAsjIndexMap = Maps.newHashMap();

        loadRefAltSjIndices(mConfig.AltSpliceJunctionSites, mRefAsjIndexMap);
    }

    @Override
    public CategoryType categoryType() { return CategoryType.ALT_SJ; }

    @Override
    public List<DataItem> extractSampleData(final String sampleId)
    {
        if(mRefAsjIndexMap.isEmpty())
            return null;

        String filename = mConfig.altSpliceJunctionFile(sampleId);

        if(!Files.exists(Paths.get(filename)))
            return null;

        List<NovelSpliceJunction> novelSpliceJunctions = NovelSpliceJunctionFile.read(filename);

        List<DataItem> dataItems = new ArrayList<>();
        Set<String> existingAsjKeys = new HashSet<>();

        for(NovelSpliceJunction novelSpliceJunction : novelSpliceJunctions)
        {
            String asjKey = NovelSpliceJunctionFile.formKey(
                    novelSpliceJunction.chromosome(), novelSpliceJunction.junctionStart(), novelSpliceJunction.junctionEnd());

            if(!mRefAsjIndexMap.containsKey(asjKey))
                continue;

            if(existingAsjKeys.contains(asjKey))
            {
                CUP_LOGGER.trace("Ignoring alt splice junction with duplicate coordinates: {}", asjKey);
                continue;
            }

            dataItems.add(new DataItem(RNA, ItemType.ALT_SJ, asjKey, novelSpliceJunction.fragmentCount()));
            existingAsjKeys.add(asjKey);
        }

        if(dataItems.isEmpty())
        {
            CUP_LOGGER.warn("sample({}) had no matching alt-SJs of the {} provided in configItem(-{})",
                    sampleId, mRefAsjIndexMap.size(), REF_ALT_SJ_SITES);
        }

        return dataItems;
    }

    protected static boolean loadRefAltSjIndices(final String filename, final Map<String,Integer> refAsjIndexMap)
    {
        try
        {
            BufferedReader fileReader = createBufferedReader(filename);

            String header = fileReader.readLine();
            String fileDelim = inferFileDelimiter(filename);
            Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(header, fileDelim);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_POS_END);

            String line = fileReader.readLine();
            int altSjIndex = 0;

            while(line != null)
            {
                final String[] items = line.split(fileDelim, -1);

                final String asjKey = NovelSpliceJunctionFile.formKey(
                        items[chrIndex], Integer.parseInt(items[posStartIndex]), Integer.parseInt(items[posEndIndex]));

                refAsjIndexMap.put(asjKey, altSjIndex++);

                line = fileReader.readLine();
            }
        }
        catch (Exception e)
        {
            CUP_LOGGER.error("Failed to read {} file({}): {}", REF_ALT_SJ_SITES, filename, e);
            System.exit(1);
        }

        return true;
    }


}
