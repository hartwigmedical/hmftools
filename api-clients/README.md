The Civic API uses cipher suites that are not enabled by default in the standard JRE. 

To retrieve Civic evidence items locally, or run the `patient-reporter` `PDFWriterTest` with `RUN_CIVIC_ANALYSIS` flag set to `true`, you need to install the unlimited strength policy JAR files. 

#### Updating cryptography policy files
 - download [Java Cryptography Extension (JCE) Unlimited Strength Jurisdiction Policy Files 8](http://www.oracle.com/technetwork/java/javase/downloads/jce8-download-2133166.html)
 - follow instructions in the included `README.txt` (**tl;dr:** copy the 2 archived jars to `JAVA_HOME/jre/lib/security`)
