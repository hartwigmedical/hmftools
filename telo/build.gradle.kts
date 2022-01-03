
plugins {
    id("com.hartwig.java-conventions")
    application
}

description = "HMF Tools - TELO"

dependencies {
    implementation(project(":hmf-common"))
    implementation(libs.tablesaw.core)
    annotationProcessor(libs.immutables.value)
    compileOnly(libs.immutables.value)
    testImplementation(libs.junit)
}

application {
    mainClass.set("com.hartwig.hmftools.telo.TeloApplication")
}
