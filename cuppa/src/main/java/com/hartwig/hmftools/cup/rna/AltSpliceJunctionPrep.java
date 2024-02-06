package com.hartwig.hmftools.cup.rna;

import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_END;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.FLD_ALT_SJ_POS_START;
import static com.hartwig.hmftools.common.rna.AltSpliceJunctionFile.formKey;
import static com.hartwig.hmftools.common.rna.RnaCommon.FLD_FRAG_COUNT;
import static com.hartwig.hmftools.common.utils.file.CommonFields.FLD_CHROMOSOME;
import static com.hartwig.hmftools.common.utils.file.FileDelimiters.inferFileDelimiter;
import static com.hartwig.hmftools.common.utils.file.FileReaderUtils.createFieldsIndexMap;
import static com.hartwig.hmftools.common.utils.file.FileWriterUtils.createBufferedReader;
import static com.hartwig.hmftools.cup.CuppaConfig.CUP_LOGGER;
import static com.hartwig.hmftools.cup.prep.DataSource.RNA;
import static com.hartwig.hmftools.cup.rna.RefAltSpliceJunctions.FLD_POS_END;
import static com.hartwig.hmftools.cup.rna.RefAltSpliceJunctions.FLD_POS_START;

import java.io.BufferedReader;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.Paths;
import java.util.List;
import java.util.Map;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.hartwig.hmftools.common.cuppa.CategoryType;
import com.hartwig.hmftools.common.rna.AltSpliceJunctionFile;
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

        List<DataItem> dataItems = Lists.newArrayList();

        final String filename = AltSpliceJunctionFile.generateFilename(mConfig.getIsofoxDataDir(sampleId), sampleId);

        if(!Files.exists(Paths.get(filename)))
            return dataItems;

        try
        {
            final List<String> lines = Files.readAllLines(Paths.get(filename));
            String fileDelim = inferFileDelimiter(filename);

            final Map<String,Integer> fieldsIndexMap = createFieldsIndexMap(lines.get(0), fileDelim);
            lines.remove(0);

            int chrIndex = fieldsIndexMap.get(FLD_CHROMOSOME);
            int posStartIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_START);
            int posEndIndex = fieldsIndexMap.get(FLD_ALT_SJ_POS_END);
            int fragCountIndex = fieldsIndexMap.get(FLD_FRAG_COUNT);

            for(String data : lines)
            {
                final String items[] = data.split(fileDelim, -1);

                String chromosome = items[chrIndex];
                int posStart = Integer.parseInt(items[posStartIndex]);
                int posEnd = Integer.parseInt(items[posEndIndex]);

                final String asjKey = formKey(chromosome, posStart, posEnd);

                Integer bucketIndex = mRefAsjIndexMap.get(asjKey);

                if(bucketIndex == null)
                    continue;

                int fragCount = Integer.parseInt(items[fragCountIndex]);

                dataItems.add(new DataItem(RNA, ItemType.ALT_SJ, asjKey, String.valueOf(fragCount)));
            }

            // CUP_LOGGER.info("loaded {} matching alt-SJs from file({})", matchedRefAltSJs, filename);
            return dataItems;

        }
        catch(IOException e)
        {
            CUP_LOGGER.error("failed to load alt splice junction file({}): {}", filename, e.toString());
            return null;
        }
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

                final String asjKey = formKey(
                        items[chrIndex], Integer.parseInt(items[posStartIndex]), Integer.parseInt(items[posEndIndex]));

                refAsjIndexMap.put(asjKey, altSjIndex++);

                line = fileReader.readLine();
            }
        }
        catch (IOException e)
        {
            CUP_LOGGER.error("failed to read RNA ref alt-SJs from file({}): {}", filename, e.toString());
            return false;
        }

        return true;
    }


}
