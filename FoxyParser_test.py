import xml.etree.ElementTree as ET

filename = input('Please enter filename: ')
# filename = 'LRG_62.xml'
print(filename)

tree = ET.parse(filename)
root = tree.getroot()

# for child in root:
#     print(child.tag, child.attrib)
#
# print('In fixed')
# for child in root[0]:
#     print(child.tag, child.attrib, child.text)
#
# print('In updateable')
# for child in root[1]:
#     print(child.tag, child.attrib)

meep = root.findall("./fixed_annotation/source/contact")[0].tag
moop = root.findall("./fixed_annotation/source/contact/email")[0].text
print(meep)
print(moop)

