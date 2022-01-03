
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - Patient Reporter"

dependencies {
    implementation(project(":hmf-common"))
    implementation(libs.log4j.slf4j.impl)
    implementation(libs.kernel)
    implementation(libs.layout)
    implementation(libs.io)
    testImplementation(libs.junit)
    testImplementation(project(path = ":hmf-common", configuration = "tests"))
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
}

application {
    mainClass.set("com.hartwig.hmftools.patientreporter.PatientReporterApplication")
}

// patient reporter tests needs 4G ram
tasks.withType<Test> {
    maxHeapSize = "4G"
}
