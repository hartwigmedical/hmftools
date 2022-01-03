plugins {
    `java-library`
}

repositories {
    mavenCentral()
}

if (!hasProperty("java.version"))
    logger.error("java.version property missing")

// set java version to the one specified in gradle.properties
val javaVersion = project.property("java.version").toString()
group = "com.hartwig"
version = "local-SNAPSHOT"
java.sourceCompatibility = JavaVersion.toVersion(javaVersion)
java.targetCompatibility = JavaVersion.toVersion(javaVersion)

tasks.withType<JavaCompile>() {
    options.encoding = "UTF-8"
}

// https://stackoverflow.com/questions/51864473/where-do-resource-files-go-in-a-gradle-project-that-builds-a-java-9-module
// the unit tests reads the resource files as normal files, so they need to be copied to the class directory
// we should fix the unit test
sourceSets {
    test {
        output.setResourcesDir(java.classesDirectory)
    }
}
