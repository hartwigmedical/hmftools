package com.hartwig.hmftools.patientdb.clinical.curators;

import static org.apache.lucene.analysis.miscellaneous.WordDelimiterGraphFilter.GENERATE_NUMBER_PARTS;
import static org.apache.lucene.analysis.miscellaneous.WordDelimiterGraphFilter.GENERATE_WORD_PARTS;
import static org.apache.lucene.analysis.miscellaneous.WordDelimiterGraphFilter.SPLIT_ON_NUMERICS;

import java.io.File;
import java.io.IOException;
import java.io.StringReader;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import com.google.common.annotations.VisibleForTesting;
import com.google.common.collect.Lists;
import com.google.common.collect.Maps;
import com.google.common.collect.Sets;
import com.hartwig.hmftools.patientdb.clinical.datamodel.CuratedDrug;
import com.hartwig.hmftools.patientdb.clinical.datamodel.ImmutableCuratedDrug;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.apache.logging.log4j.util.Strings;
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
import org.apache.lucene.search.highlight.QueryTermExtractor;
import org.apache.lucene.search.highlight.WeightedTerm;
import org.apache.lucene.search.spell.Dictionary;
import org.apache.lucene.search.spell.HighFrequencyDictionary;
import org.apache.lucene.search.spell.SpellChecker;
import org.apache.lucene.store.Directory;
import org.apache.lucene.store.RAMDirectory;
import org.jetbrains.annotations.NotNull;

public class TreatmentCurator implements CleanableCurator {

    private static final Logger LOGGER = LogManager.getLogger(TreatmentCurator.class);
    private static final String FIELD_DELIMITER = "\t";
    private static final String SYNONYM_DELIMITER = ",";

    private static final String CANONICAL_DRUG_NAME_FIELD = "canonicalDrugName";
    private static final String DRUG_NAME_FIELD = "drugName";
    private static final String DRUG_TYPE_FIELD = "drugType";
    private static final String DRUG_MECHANISM_FIELD = "drugMechanism";
    private static final String DRUG_SYNONYMS_FIELD = "drugSynonyms";

    private static final int NUM_HITS = 20;
    private static final int MAX_SHINGLES = 10;
    private static final float SPELLCHECK_ACCURACY = .85f;

    @NotNull
    private final Map<String, String> unusedTokenizedTermToEntryMap;
    @NotNull
    private final SpellChecker spellChecker;
    @NotNull
    private final IndexSearcher indexSearcher;

    @VisibleForTesting
    public TreatmentCurator(@NotNull String treatmentMappingTsv) throws IOException {
        List<DrugEntry> drugEntries = readEntries(treatmentMappingTsv);
        Directory index = createIndex(drugEntries);
        IndexReader reader = DirectoryReader.open(index);

        spellChecker = createIndexSpellchecker(index);
        indexSearcher = new IndexSearcher(reader);
        unusedTokenizedTermToEntryMap = extractUnusedTokenizedTermToEntryMap(drugEntries);
    }

    @NotNull
    private static Map<String, String> extractUnusedTokenizedTermToEntryMap(@NotNull Iterable<DrugEntry> drugEntries) throws IOException {
        Map<String, String> uniqueTokenizedTermToEntryMap = Maps.newHashMap();

        for (DrugEntry drugEntry : drugEntries) {
            for (String synonym : drugEntry.synonyms()) {
                if (uniqueTokenizedTermToEntryMap.containsValue(synonym)) {
                    LOGGER.warn("Drug synonym already included in search terms: {}", synonym);
                } else {
                    uniqueTokenizedTermToEntryMap.put(toTokenizedString(synonym), synonym);
                }
            }
            if (drugEntry.synonyms().isEmpty()) {
                if (uniqueTokenizedTermToEntryMap.containsValue(drugEntry.canonicalName())) {
                    LOGGER.warn("Drug canonical name already included in search terms: {}", drugEntry.canonicalName());
                } else {
                    uniqueTokenizedTermToEntryMap.put(toTokenizedString(drugEntry.canonicalName()), drugEntry.canonicalName());
                }
            }
        }

        return uniqueTokenizedTermToEntryMap;
    }

    @NotNull
    private static String toTokenizedString(@NotNull String searchTerm) throws IOException {
        List<String> result = Lists.newArrayList();
        TokenStream stream = indexAnalyzer().tokenStream(DRUG_NAME_FIELD, new StringReader(searchTerm));
        stream.reset();
        while (stream.incrementToken()) {
            result.add(stream.getAttribute(CharTermAttribute.class).toString());
        }
        assert result.size() == 1;
        return result.get(0);
    }

    @NotNull
    @Override
    public Set<String> unusedSearchTerms() {
        return Sets.newHashSet(unusedTokenizedTermToEntryMap.values());
    }

    @NotNull
    private static List<DrugEntry> readEntries(@NotNull String treatmentMappingTsv) throws IOException {
        List<DrugEntry> drugEntries = Lists.newArrayList();

        List<String> lines = Files.readAllLines(new File(treatmentMappingTsv).toPath());

        // Skip header
        for (String line : lines.subList(1, lines.size())) {
            String[] parts = line.split(FIELD_DELIMITER);
            String drugCanonicalName = parts[0];
            String treatmentType = parts[1];
            String treatmentMechanism = parts.length > 2 ? parts[2] : Strings.EMPTY;
            String synonymsField = parts.length > 3 ? parts[3] :Strings.EMPTY;

            drugEntries.add(ImmutableDrugEntry.builder()
                    .canonicalName(drugCanonicalName)
                    .type(treatmentType)
                    .treatmentMechanism(treatmentMechanism)
                    .synonyms(toSynonyms(synonymsField))
                    .build());
        }
        return drugEntries;
    }

    @NotNull
    @VisibleForTesting
    static List<String> toSynonyms(@NotNull String synonymField) {
        if (synonymField.isEmpty()) {
            return Lists.newArrayList();
        }

        List<String> synonyms = Lists.newArrayList();
        // Split on comma but ignore comma's within quotes
        for (String token : synonymField.split(",(?=(?:[^\"]*\"[^\"]*\")*[^\"]*$)", -1)) {
            String trimmed = token.trim();
            if (trimmed.startsWith("\"") && trimmed.endsWith("\"")) {
                synonyms.add(trimmed.substring(1, trimmed.length() - 1));
            } else {
                synonyms.add(trimmed);
            }
        }
        return synonyms;
    }

    @NotNull
    public List<CuratedDrug> search(@NotNull String searchTerm) {
        Optional<CuratedDrug> matchedDrug = matchSingle(searchTerm);
        if (!matchedDrug.isPresent()) {
            return matchMultiple(searchTerm);
        } else {
            return Lists.newArrayList(matchedDrug.get());
        }
    }

    @NotNull
    Optional<CuratedDrug> matchSingle(@NotNull String searchTerm) {
        Analyzer analyzer = spellcheckAnalyzer(spellChecker);
        Query query = new QueryParser(DRUG_NAME_FIELD, analyzer).createPhraseQuery(DRUG_NAME_FIELD, searchTerm);
        try {
            ScoreDoc[] hits = indexSearcher.search(query, NUM_HITS).scoreDocs;

            for (WeightedTerm term : QueryTermExtractor.getTerms(query)) {
                unusedTokenizedTermToEntryMap.remove(term.getTerm());
            }

            if (hits.length == 1) {
                Document searchResult = indexSearcher.doc(hits[0].doc);
                return Optional.of(ImmutableCuratedDrug.of(searchResult.get(CANONICAL_DRUG_NAME_FIELD),
                        searchResult.get(DRUG_TYPE_FIELD),
                        searchResult.get(DRUG_MECHANISM_FIELD),
                        searchTerm));
            }

            return Optional.empty();
        } catch (IOException | NullPointerException exception) {
            // NullPointerException raised when query is null
            // (happens when search term contains only separator characters like whitespace, commas, punctuation, etc)
            LOGGER.error("Caught {} in treatment curation for search term {}. Error message: {} ",
                    exception.getClass().getName(),
                    searchTerm,
                    exception.getMessage());
            return Optional.empty();
        }
    }

    @NotNull
    List<CuratedDrug> matchMultiple(@NotNull String searchTerm) {
        Map<SearchToken, CuratedDrug> matchedTokens = Maps.newHashMap();
        List<SearchToken> searchTokens = generateSearchTokens(searchTerm);
        for (SearchToken searchToken : searchTokens) {
            if (matchedTokens.keySet()
                    .stream()
                    .noneMatch(token -> (token.startOffset() <= searchToken.startOffset() && searchToken.startOffset() <= token.endOffset())
                            || (token.startOffset() <= searchToken.endOffset() && searchToken.endOffset() <= token.endOffset()))) {
                Optional<CuratedDrug> matchedTreatment = matchSingle(searchToken.term());
                matchedTreatment.ifPresent(curatedTreatment -> matchedTokens.put(searchToken, curatedTreatment));
            }
        }
        return Lists.newArrayList(matchedTokens.values());
    }

    @NotNull
    private static List<SearchToken> generateSearchTokens(@NotNull String searchTerm) {
        Set<SearchToken> searchTokens = Sets.newHashSet();
        TokenStream tokenStream = getSpellCheckedShingleStream(searchTerm);
        try {
            tokenStream.reset();

            while (tokenStream.incrementToken()) {
                String searchToken = tokenStream.getAttribute(CharTermAttribute.class).toString();
                OffsetAttribute offsetAttribute = tokenStream.getAttribute(OffsetAttribute.class);
                searchTokens.add(ImmutableSearchToken.of(searchToken, offsetAttribute.startOffset(), offsetAttribute.endOffset()));
            }
            tokenStream.end();
            tokenStream.close();
            return Lists.newArrayList(searchTokens);
        } catch (IOException exception) {
            LOGGER.warn("Caught IOException in treatment curation: {}", exception.getMessage());
            return Lists.newArrayList();
        }
    }

    @NotNull
    private static TokenStream getSpellCheckedShingleStream(@NotNull String searchTerm) {
        StringReader reader = new StringReader(searchTerm);
        Analyzer analyzer = createShingleAnalyzer(MAX_SHINGLES);
        return analyzer.tokenStream(DRUG_NAME_FIELD, reader);
    }

    @NotNull
    private static Directory createIndex(@NotNull List<DrugEntry> drugEntries) throws IOException {
        Directory drugIndex = new RAMDirectory();
        IndexWriter indexWriter = createIndexWriter(drugIndex);
        for (DrugEntry drugEntry : drugEntries) {
            indexDrugEntry(indexWriter, drugEntry);
        }
        indexWriter.close();
        return drugIndex;
    }

    @NotNull
    private static IndexWriter createIndexWriter(@NotNull Directory directory) throws IOException {
        Analyzer analyzer = indexAnalyzer();
        IndexWriterConfig config = new IndexWriterConfig(analyzer);
        return new IndexWriter(directory, config);
    }

    private static void indexDrugEntry(@NotNull IndexWriter writer, @NotNull DrugEntry drugEntry) throws IOException {
        Document document = new Document();
        drugEntry.synonyms().forEach(synonym -> {
            document.add(new TextField(DRUG_NAME_FIELD, synonym, Field.Store.NO));
            document.add(new TextField(DRUG_SYNONYMS_FIELD, synonym, Field.Store.YES));
        });
        document.add(new TextField(DRUG_NAME_FIELD, drugEntry.canonicalName(), Field.Store.NO));
        document.add(new TextField(DRUG_SYNONYMS_FIELD, drugEntry.canonicalName(), Field.Store.YES));
        document.add(new StringField(DRUG_TYPE_FIELD, drugEntry.type(), Field.Store.YES));
        document.add(new TextField(DRUG_MECHANISM_FIELD, drugEntry.treatmentMechanism(), Field.Store.YES));
        document.add(new TextField(CANONICAL_DRUG_NAME_FIELD, drugEntry.canonicalName(), Field.Store.YES));
        writer.addDocument(document);
    }

    @NotNull
    private static SpellChecker createIndexSpellchecker(@NotNull Directory index) throws IOException {
        Directory spellCheckerDirectory = new RAMDirectory();
        IndexReader indexReader = DirectoryReader.open(index);
        Analyzer analyzer = new SimpleAnalyzer();
        IndexWriterConfig config = new IndexWriterConfig(analyzer);
        Dictionary dictionary = new HighFrequencyDictionary(indexReader, DRUG_SYNONYMS_FIELD, 0.0f);
        SpellChecker spellChecker = new SpellChecker(spellCheckerDirectory);

        spellChecker.indexDictionary(dictionary, config, false);
        spellChecker.setAccuracy(SPELLCHECK_ACCURACY);
        return spellChecker;
    }

    @NotNull
    private static Analyzer createShingleAnalyzer(int maxShingles) {
        return new Analyzer() {
            @Override
            protected TokenStreamComponents createComponents(@NotNull String field) {
                Tokenizer source = new WhitespaceTokenizer();
                source.setReader(new StringReader(field));
                ShingleFilter shingleFilter = new ShingleFilter(defaultTokenFilter(source), maxShingles);
                shingleFilter.setOutputUnigrams(true);
                return new TokenStreamComponents(source, shingleFilter);
            }
        };
    }

    @NotNull
    private static Analyzer spellcheckAnalyzer(@NotNull SpellChecker spellChecker) {
        return new Analyzer() {
            @Override
            protected TokenStreamComponents createComponents(@NotNull String field) {
                Tokenizer source = new WhitespaceTokenizer();
                source.setReader(new StringReader(field));
                SpellCheckerTokenFilter spellCheckFilter = new SpellCheckerTokenFilter(defaultTokenFilter(source), spellChecker);
                TokenFilter concatenatingFilter = new ConcatenatingFilter(spellCheckFilter, ' ');
                return new TokenStreamComponents(source, concatenatingFilter);
            }
        };
    }

    @NotNull
    private static Analyzer wordDelimiterAnalyzer() {
        return new Analyzer() {
            @Override
            protected TokenStreamComponents createComponents(@NotNull String field) {
                Tokenizer source = new WhitespaceTokenizer();
                source.setReader(new StringReader(field));
                return new TokenStreamComponents(source, defaultTokenFilter(source));
            }
        };
    }

    @NotNull
    private static Analyzer concatenatingAnalyzer() {
        return new Analyzer() {
            @Override
            protected TokenStreamComponents createComponents(@NotNull String field) {
                Tokenizer source = new WhitespaceTokenizer();
                source.setReader(new StringReader(field));
                TokenFilter concatenatingFilter = new ConcatenatingFilter(defaultTokenFilter(source), ' ');
                return new TokenStreamComponents(source, concatenatingFilter);
            }
        };
    }

    @NotNull
    private static TokenFilter defaultTokenFilter(@NotNull Tokenizer source) {
        TokenFilter filteredSource = new LowerCaseFilter(source);
        return new WordDelimiterGraphFilter(filteredSource, SPLIT_ON_NUMERICS | GENERATE_WORD_PARTS | GENERATE_NUMBER_PARTS, null);
    }

    @NotNull
    private static Analyzer indexAnalyzer() {
        Map<String, Analyzer> fieldAnalyzers = Maps.newHashMap();
        fieldAnalyzers.put(DRUG_NAME_FIELD, concatenatingAnalyzer());
        return new PerFieldAnalyzerWrapper(wordDelimiterAnalyzer(), fieldAnalyzers);
    }
}
