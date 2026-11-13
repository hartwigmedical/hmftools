package com.hartwig.hmftools.compar;

import static com.hartwig.hmftools.common.utils.config.CommonConfig.SAMPLE;
import static com.hartwig.hmftools.compar.ComparConfig.CATEGORIES;
import static com.hartwig.hmftools.compar.ComparConfig.MATCH_LEVEL;
import static com.hartwig.hmftools.compar.common.CategoryType.ALL_CATEGORIES;

import static org.junit.Assert.assertEquals;
import static org.junit.Assert.assertNotNull;
import static org.junit.Assert.assertNull;

import java.io.File;
import java.io.IOException;
import java.util.List;

import com.hartwig.hmftools.common.utils.config.ConfigBuilder;
import com.hartwig.hmftools.compar.common.MatchLevel;
import com.hartwig.hmftools.compar.common.field.FieldCheck;
import com.hartwig.hmftools.compar.common.field.FieldInfo;

import org.junit.Rule;
import org.junit.Test;
import org.junit.rules.TemporaryFolder;

public class FieldConfigFileTest
{
    @Rule
    public final TemporaryFolder temp = new TemporaryFolder();

    @Test
    public void writtenFieldConfigLoadsBackUnchanged() throws IOException
    {
        for(MatchLevel matchLevel : MatchLevel.values())
        {
            ConfigBuilder configBuilder = new ConfigBuilder();
            ComparConfig.addConfig(configBuilder);
            configBuilder.setValue(CATEGORIES, ALL_CATEGORIES);
            configBuilder.setValue(MATCH_LEVEL, matchLevel.toString());
            configBuilder.setValue(SAMPLE, "TEST");

            List<ItemComparer> comparers = ComparerUtils.buildComparers(new ComparConfig(configBuilder), new FieldCheckCache());

            File fieldConfigFile = temp.newFile();
            FieldConfigFile.write(fieldConfigFile.getAbsolutePath(), comparers);

            FieldCheckCache fieldCheckCache = new FieldCheckCache();
            fieldCheckCache.loadOverrides(fieldConfigFile.getAbsolutePath());

            for(ItemComparer comparer : comparers)
            {
                for(FieldInfo field : comparer.fieldsList())
                {
                    if(field.DisplayOnly)
                        continue;

                    FieldCheck loaded = fieldCheckCache.getCategoryFieldCheckOverrides(comparer.category()).get(field.Name);
                    String context = matchLevel + " " + comparer.category() + " " + field.Name;

                    assertNotNull(context, loaded);
                    assertEquals(context, field.FieldCheck.IsCompared, loaded.IsCompared);
                    assertThreshold(context, field.FieldCheck.AbsDiffThreshold, loaded.AbsDiffThreshold);
                    assertThreshold(context, field.FieldCheck.PercDiffThreshold, loaded.PercDiffThreshold);
                }
            }
        }
    }

    private static void assertThreshold(final String context, final Double expected, final Double actual)
    {
        if(expected == null)
            assertNull(context, actual);
        else
            assertEquals(context, expected, actual, 1e-9);
    }
}
