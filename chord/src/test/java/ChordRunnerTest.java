import java.io.File;
import java.io.IOException;
import java.util.List;

import com.google.common.io.Resources;
import com.hartwig.hmftools.chord.ChordConfig;
import com.hartwig.hmftools.chord.ChordRunner;

import org.apache.commons.io.FileUtils;
import org.junit.After;
import org.junit.Ignore;
import org.junit.Test;

@Ignore
public class ChordRunnerTest
{
    private static final String OUTPUT_DIR = System.getProperty("java.io.tmpdir") + "/chord_output/";

    @After
    public void teardown() throws IOException
    {
        FileUtils.deleteDirectory(new File(OUTPUT_DIR));
    }

    @Test
    public void canRunChordOnEmptyVcfs()
    {
        String sampleId = "EMPTY_SAMPLE";
        String inputDir = Resources.getResource("vcf/").getPath();
        String chordToolDir = "/Users/lnguyen/Hartwig/hartwigmedical/hmftools/chord/src/main/R/";

        new File(OUTPUT_DIR).mkdirs();

        ChordConfig config = new ChordConfig.Builder()
                .sampleIds(List.of(sampleId))
                .purpleDir(inputDir)
                .outputDir(OUTPUT_DIR)
                .chordToolDir(chordToolDir)
                .build();

        ChordRunner runner = new ChordRunner(config);
        runner.run();
    }
}
