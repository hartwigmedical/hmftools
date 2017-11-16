package com.hartwig.hmftools.patientdb.curators;

import static org.apache.lucene.analysis.miscellaneous.WordDelimiterGraphFilter.GENERATE_NUMBER_PARTS;
import static org.apache.lucene.analysis.miscellaneous.WordDelimiterGraphFilter.GENERATE_WORD_PARTS;
import static org.apache.lucene.analysis.miscellaneous.WordDelimiterGraphFilter.SPLIT_ON_NUMERICS;

import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.nio.charset.Charset;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientdb.data.CuratedTreatment;
import com.hartwig.hmftools.patientdb.data.ImmutableCuratedTreatment;

import org.apache.commons.csv.CSVFormat;
import org.apache.commons.csv.CSVParser;
import org.apache.commons.csv.CSVRecord;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.lucene.analysis.Analyzer;
import org.apache.lucene.analysis.LowerCaseFilter;
import org.apache.lucene.analysis.TokenFilter;
import org.apache.lucene.analysis.TokenStream;
import org.apache.lucene.analysis.Tokenizer;
import org.apache.lucene.analysis.core.SimpleAnalyzer;
import org.apache.lucene.analysis.core.WhitespaceTokenizer;
import org.apache.lucene.analysis.miscellaneous.PerFieldAnalyzerWrapper;
import org.apache.lucene.analysis.miscellaneous.WordDelimiterGraphFilter;
import org.apache.lucene.analysis.shingle.ShingleFilter;
import org.apache.lucene.analysis.tokenattributes.CharTermAttribute;
import org.apache.lucene.analysis.tokenattributes.OffsetAttribute;
import org.apache.lucene.document.Document;
import org.apache.lucene.document.Field;
import org.apache.lucene.document.StringField;
import org.apache.lucene.document.TextField;
import org.apache.lucene.index.DirectoryReader;
import org.apache.lucene.index.IndexReader;
import org.apache.lucene.index.IndexWriter;
import org.apache.lucene.index.IndexWriterConfig;
import org.apache.lucene.queryparser.classic.QueryParser;
import org.apache.lucene.search.IndexSearcher;
import org.apache.lucene.search.Query;
import org.apache.lucene.search.ScoreDoc;
import org.apache.lucene.search.spell.Dictionary;
import org.apache.lucene.search.spell.HighFrequencyDictionary;
import org.apache.lucene.search.spell.SpellChecker;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.RAMDirectory;
import org.jetbrains.annotations.NotNull;

public class TreatmentCurator {
    private static final Logger LOGGER = LogManager.getLogger(TreatmentCurator.class);
    private static final String DRUG_TERMS_FIELD = "drugTerms";
    private static final String DRUG_NAME_FIELD = "drugName";
    private static final String CANONICAL_DRUG_NAME_FIELD = "canonicalDrugName";
    private static final String DRUG_TYPE_FIELD = "drugType";
    private static final String CANONICAL_DRUG_NAME_CSV_FIELD = "drug";
    private static final String DRUG_TYPE_CSV_FIELD = "type";
    private static final String OTHER_DRUG_NAMES_CSV_FIELD = "other_names";
    private static final int NUM_HITS = 20;
    private static final int MAX_SHINGLES = 10;
    private static final float SPELLCHECK_ACCURACY = .85f;
    private final SpellChecker spellChecker;
    private final IndexSearcher indexSearcher;

    public TreatmentCurator(@NotNull final String mappingCsv) throws IOException {
        final Directory index = createIndex(mappingCsv);
        final IndexReader reader = DirectoryReader.open(index);
        spellChecker = createIndexSpellchecker(index);
        indexSearcher = new IndexSearcher(reader);
    }

    @NotNull
    public List<CuratedTreatment> search(@NotNull final String searchTerm) throws IOException {
        final Optional<CuratedTreatment> matchedTreatment = matchSingle(searchTerm);
        if (!matchedTreatment.isPresent()) {
            LOGGER.warn("Failed to match search term: {}. attempting to match multiple treatments", searchTerm);
            final List<CuratedTreatment> matchedTreatments = matchMultiple(searchTerm);
            if (!matchedTreatments.isEmpty()) {
                LOGGER.info("Matched multiple treatments {} to {}", searchTerm,
                        matchedTreatments.stream().map(CuratedTreatment::name).collect(Collectors.toList()));
            }
            return matchedTreatments;
        } else {
            return Lists.newArrayList(matchedTreatment.get());
        }
    }

    @NotNull
    Optional<CuratedTreatment> matchSingle(@NotNull final String searchTerm) throws IOException {
        final Analyzer analyzer = spellcheckAnalyzer(spellChecker);
        final Query query = new QueryParser(DRUG_NAME_FIELD, analyzer).createPhraseQuery(DRUG_NAME_FIELD, searchTerm);
        final ScoreDoc[] hits = indexSearcher.search(query, NUM_HITS).scoreDocs;
        if (hits.length == 1) {
            final Document searchResult = indexSearcher.doc(hits[0].doc);
            return Optional.of(ImmutableCuratedTreatment.of(searchResult.get(CANONICAL_DRUG_NAME_FIELD), searchResult.get(DRUG_TYPE_FIELD),
                    searchTerm));
        }
        return Optional.empty();
    }

    @NotNull
    List<CuratedTreatment> matchMultiple(@NotNull final String searchTerm) throws IOException {
        final HashMap<SearchToken, CuratedTreatment> matchedTokens = Maps.newHashMap();
        final List<SearchToken> searchTokens = generateSearchTokens(searchTerm);
        for (final SearchToken searchToken : searchTokens) {
            if (matchedTokens.keySet()
                    .stream()
                    .noneMatch(token -> (token.startOffset() <= searchToken.startOffset() && searchToken.startOffset() <= token.endOffset())
                            || (token.startOffset() <= searchToken.endOffset() && searchToken.endOffset() <= token.endOffset()))) {
                final Optional<CuratedTreatment> matchedTreatment = matchSingle(searchToken.term());
                matchedTreatment.ifPresent(curatedTreatment -> matchedTokens.put(searchToken, curatedTreatment));
            }
        }
        return Lists.newArrayList(matchedTokens.values());
    }

    @NotNull
    private static List<SearchToken> generateSearchTokens(@NotNull final String searchTerm) throws IOException {
        final Set<SearchToken> searchTokens = Sets.newHashSet();
        final TokenStream tokenStream = getSpellCheckedShingleStream(searchTerm);
        tokenStream.reset();
        while (tokenStream.incrementToken()) {
            final String searchToken = tokenStream.getAttribute(CharTermAttribute.class).toString();
            final OffsetAttribute offsetAttribute = tokenStream.getAttribute(OffsetAttribute.class);
            searchTokens.add(ImmutableSearchToken.of(searchToken, offsetAttribute.startOffset(), offsetAttribute.endOffset()));
        }
        tokenStream.end();
        tokenStream.close();
        return searchTokens.stream().sorted(Comparator.comparing(SearchToken::length).reversed().thenComparing(SearchToken::startOffset)).collect(Collectors.toList());
    }

    @NotNull
    private static TokenStream getSpellCheckedShingleStream(@NotNull final String searchTerm) {
        StringReader reader = new StringReader(searchTerm);
        final Analyzer analyzer = createShingleAnalyzer(MAX_SHINGLES);
        return analyzer.tokenStream(DRUG_NAME_FIELD, reader);
    }

    @NotNull
    private static Directory createIndex(@NotNull final String mappingCsv) throws IOException {
        final Directory treatmentIndex = new RAMDirectory();
        final CSVParser parser = CSVParser.parse(new File(mappingCsv), Charset.defaultCharset(), CSVFormat.DEFAULT.withHeader());
        final IndexWriter indexWriter = createIndexWriter(treatmentIndex);
        for (final CSVRecord record : parser) {
            indexRecord(indexWriter, record);
        }
        indexWriter.close();
        return treatmentIndex;
    }

    @NotNull
    private static IndexWriter createIndexWriter(@NotNull final Directory directory) throws IOException {
        final Analyzer analyzer = indexAnalyzer();
        final IndexWriterConfig config = new IndexWriterConfig(analyzer);
        return new IndexWriter(directory, config);
    }

    @NotNull
    private static SpellChecker createIndexSpellchecker(@NotNull final Directory index) throws IOException {
        final Directory spellCheckerDirectory = new RAMDirectory();
        final IndexReader indexReader = DirectoryReader.open(index);
        final Analyzer analyzer = new SimpleAnalyzer();
        final IndexWriterConfig config = new IndexWriterConfig(analyzer);
        final Dictionary dictionary = new HighFrequencyDictionary(indexReader, DRUG_TERMS_FIELD, 0.0f);
        final SpellChecker spellChecker = new SpellChecker(spellCheckerDirectory);

        spellChecker.indexDictionary(dictionary, config, false);
        spellChecker.setAccuracy(SPELLCHECK_ACCURACY);
        return spellChecker;
    }

    private static void indexRecord(@NotNull final IndexWriter writer, @NotNull final CSVRecord record) throws IOException {
        final String canonicalName = record.get(CANONICAL_DRUG_NAME_CSV_FIELD);
        final String drugType = record.get(DRUG_TYPE_CSV_FIELD);
        final String otherNamesString = record.get(OTHER_DRUG_NAMES_CSV_FIELD);
        final List<String> drugNames = Lists.newArrayList();
        if (!otherNamesString.isEmpty()) {
            final CSVParser otherNamesParser = CSVParser.parse(otherNamesString, CSVFormat.DEFAULT);
            for (final CSVRecord otherNames : otherNamesParser) {
                for (final String name : otherNames) {
                    drugNames.add(name.trim());
                }
            }
        }
        indexEntry(writer, drugNames, drugType, canonicalName.trim());
    }

    private static void indexEntry(@NotNull final IndexWriter writer, @NotNull final List<String> names, @NotNull final String type,
            @NotNull final String canonicalName) throws IOException {
        final Document document = new Document();
        names.forEach(name -> {
            document.add(new TextField(DRUG_NAME_FIELD, name, Field.Store.NO));
            document.add(new TextField(DRUG_TERMS_FIELD, name, Field.Store.YES));
        });
        document.add(new TextField(DRUG_NAME_FIELD, canonicalName, Field.Store.NO));
        document.add(new TextField(DRUG_TERMS_FIELD, canonicalName, Field.Store.YES));
        document.add(new StringField(DRUG_TYPE_FIELD, type, Field.Store.YES));
        document.add(new TextField(CANONICAL_DRUG_NAME_FIELD, canonicalName, Field.Store.YES));
        writer.addDocument(document);
    }

    @NotNull
    private static Analyzer createShingleAnalyzer(final int maxShingles) {
        return new Analyzer() {
            @Override
            protected TokenStreamComponents createComponents(@NotNull final String field) {
                final Tokenizer source = new WhitespaceTokenizer();
                source.setReader(new StringReader(field));
                final ShingleFilter shingleFilter = new ShingleFilter(defaultTokenFilter(source), maxShingles);
                shingleFilter.setOutputUnigrams(true);
                return new TokenStreamComponents(source, shingleFilter);
            }
        };
    }

    @NotNull
    private static Analyzer spellcheckAnalyzer(@NotNull final SpellChecker spellChecker) {
        return new Analyzer() {
            @Override
            protected TokenStreamComponents createComponents(@NotNull final String field) {
                final Tokenizer source = new WhitespaceTokenizer();
                source.setReader(new StringReader(field));
                final SpellCheckerTokenFilter spellCheckFilter = new SpellCheckerTokenFilter(defaultTokenFilter(source), spellChecker);
                final TokenFilter concatenatingFilter = new ConcatenatingFilter(spellCheckFilter, ' ');
                return new TokenStreamComponents(source, concatenatingFilter);
            }
        };
    }

    @NotNull
    private static Analyzer wordDelimiterAnalyzer() {
        return new Analyzer() {
            @Override
            protected TokenStreamComponents createComponents(@NotNull final String field) {
                final Tokenizer source = new WhitespaceTokenizer();
                source.setReader(new StringReader(field));
                return new TokenStreamComponents(source, defaultTokenFilter(source));
            }
        };
    }

    @NotNull
    private static Analyzer concatenatingAnalyzer() {
        return new Analyzer() {
            @Override
            protected TokenStreamComponents createComponents(@NotNull final String field) {
                final Tokenizer source = new WhitespaceTokenizer();
                source.setReader(new StringReader(field));
                final TokenFilter concatenatingFilter = new ConcatenatingFilter(defaultTokenFilter(source), ' ');
                return new TokenStreamComponents(source, concatenatingFilter);
            }
        };
    }

    @NotNull
    private static TokenFilter defaultTokenFilter(@NotNull final Tokenizer source) {
        final TokenFilter filteredSource = new LowerCaseFilter(source);
        return new WordDelimiterGraphFilter(filteredSource, SPLIT_ON_NUMERICS | GENERATE_WORD_PARTS | GENERATE_NUMBER_PARTS, null);
    }

    @NotNull
    private static Analyzer indexAnalyzer() {
        final Map<String, Analyzer> fieldAnalyzers = Maps.newHashMap();
        fieldAnalyzers.put(DRUG_NAME_FIELD, concatenatingAnalyzer());
        return new PerFieldAnalyzerWrapper(wordDelimiterAnalyzer(), fieldAnalyzers);
    }
}
