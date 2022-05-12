Repository for placing the course manuals - these should be markdown or PDF documents

Converting between markdown to PDF can be perfored using pandoc, for example:

pandoc --pdf-engine=xelatex module_artemis.md -o module_artemis.pdf

sd21 - also can do it with the convert tool in atom - this seems to work a bit better as pandoc was misformatting some documents based on strange fonts
