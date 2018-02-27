package com.hartwig.hmftools.bachelor.predicates;

import static org.junit.Assert.assertFalse;
import static org.junit.Assert.assertTrue;

import java.io.IOException;
import java.nio.file.Paths;
import java.util.function.Predicate;

import javax.xml.bind.JAXBException;

import com.google.common.collect.Sets;
import com.google.common.io.Resources;
import com.hartwig.hmftools.bachelor.BachelorSchema;
import com.hartwig.hmftools.bachelor.VariantModel;

import nl.hartwigmedicalfoundation.bachelor.Program;

import org.jetbrains.annotations.NotNull;
import org.junit.Before;
import org.junit.Test;
import org.xml.sax.SAXException;

import htsjdk.variant.variantcontext.VariantContext;
import htsjdk.variant.vcf.VCFCodec;
import htsjdk.variant.vcf.VCFHeader;
import htsjdk.variant.vcf.VCFHeaderVersion;

public class BlacklistPredicateTest {

    private final static String BLACKLIST_XML = Resources.getResource("blacklist.xml").getPath();
    private VCFCodec codec;

    private static final String START = "13\t32959149\t.\tC\tT\t.\tPASS\tAC=0;AF=0;AN=0;NT=ref;QSS=40;QSS_NT=40;SGT=CC->CT;TQSS=2;";
    private static final String ANN =
            "ANN=T|intron_variant|MODIFIER|BRCA2|ENSG00000139618|transcript|ENST00000380152|protein_coding|24/26|c.9256+4867C>A|p.Ser430Glx|1307/3444|1288/2043|430/680;";
    private static final String END = "MAPPABILITY=1.000000\tGT:AD:DP\t0/1:98,21:121";
    private static final String SAMPLE = "sample";

    private Predicate<VariantModel> blacklistPredicate;

    @Before
    public void setup() throws JAXBException, SAXException {
        BachelorSchema schema = BachelorSchema.make();
        final Program program = schema.processXML(Paths.get(BLACKLIST_XML));
        blacklistPredicate = new BlacklistPredicate(Sets.newHashSet("ENST00000380152"), program.getBlacklist());

        VCFHeader header = new VCFHeader(Sets.newHashSet(), Sets.newHashSet(SAMPLE));
        codec = new VCFCodec();
        codec.setVCFHeader(header, VCFHeaderVersion.VCF4_2);
    }

    @Test
    public void testNotInBlacklist() throws JAXBException, IOException, SAXException {
        VariantModel model = createModel(START + ANN + END);
        assertFalse(blacklistPredicate.test(model));
    }

    @Test
    public void testHGVSc() throws JAXBException, IOException, SAXException {
        VariantModel model = createModel(START + ANN.replace("c.9256+4867C>A", "c.9256+4867C>T") + END);
        assertTrue(blacklistPredicate.test(model));
    }

    @Test
    public void testHGVSp() throws JAXBException, IOException, SAXException {
        VariantModel model = createModel(START + ANN.replace("p.Ser430Glx", "p.Ser430Gly") + END);
        assertTrue(blacklistPredicate.test(model));
    }

    @Test
    public void testMinCodon() throws JAXBException, IOException, SAXException {
        VariantModel model = createModel(START + ANN.replace("430/680", "1310/1680") + END);
        assertTrue(blacklistPredicate.test(model));
    }

    @NotNull
    private VariantModel createModel(@NotNull final String line) {
        VariantContext context = codec.decode(line);
        return new VariantModel(SAMPLE, context);
    }

}
