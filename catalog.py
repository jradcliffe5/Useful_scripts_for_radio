#!/usr/bin/env python

# Display the contents of AIPS catalog files.
# 
# Stephen Bourke
# Version: 20090203
# 2006 -

"""Print the contents of AIPS data directories (AIPS disks.)"""

import commands, dircache, glob, optparse, os, pwd, re, string, struct, sys 

def int2date(date):
	"""Convert the AIPS catalog date value to a string."""
	# Date is of the form 0xYYMMDD
	months = ['ERR', 'JAN', 'FEB', 'MAR', 'APR', 'MAY', 'JUN', \
			 'JUL', 'AUG', 'SEP', 'OCT', 'NOV', 'DEC']
	dd = date & 0xff
	mm = (date & 0xff00) >> 8
	yy = ((date & 0xff0000) >> 16) + 1900
	return "%02d-%s-%d" % (dd, months[mm], yy)

def int2time(time):
	"""Convert the AIPS catalog time value to a string."""
	# Time is of the form 0xHHMMSS
	ss = time & 0xff
	mm = (time & 0xff00) >> 8
	hh = (time & 0xff0000) >> 16
	return "%02d:%02d:%02d" % (hh, mm, ss)

def int2stat(stat):
	"""Convert the AIPS catalog status value to a string."""
	# See Going AIPS 5.4.4
	if stat > 0:
		return 'READ'
	elif stat < 0:
		return 'WRIT'
	else:
		return '    '

def parse_catalog(catalog_file):
	"""Return a list of catalog entry dictionaries."""
	# See Going AIPS 5.4
	catalog_format = 'iiiii12s6s2s' # Binary format of entries in catalog
	block_size = 1024
	record_size = struct.calcsize(catalog_format)
	catalog_entries = []
	
	catalog = open(catalog_file, 'rb')
	# Seek past the Header block
	catalog.seek(block_size)
	entry_raw = catalog.read(record_size)
	
	slot = 1
	while len(entry_raw) == record_size:
		id, stat, date, time, seq, inname, inclass, type = \
			struct.unpack(catalog_format, entry_raw)
		if id != -1:
			catalog_entries.append({'slot': slot, 'id': id, \
						'inname': inname, \
						'inclass': inclass, \
						'seq': seq, 'type': type, \
						'date': int2date(date), \
						'time': int2time(time), \
						'stat': int2stat(stat)})
		# Catalog entries are block aligned
		position_in_block = catalog.tell() % block_size
		if position_in_block + record_size > block_size:
			catalog.seek(block_size - position_in_block, 1)
		entry_raw = catalog.read(record_size)
		slot += 1
	return catalog_entries

def entry2str(entry):
	"""Convert a catalog entry to a PCAT-like string."""
	return '%5d%5d %s.%s.%5d %s %s %s %s' % \
			(entry['slot'], entry['id'], entry['inname'], \
			 entry['inclass'], entry['seq'], entry['type'], \
			 entry['date'], entry['time'], entry['stat'])

def to_basex(num, chars, width=0):
	"""Convert a number into arbitary base."""
	base = len(chars)
	
	hex_str = ''
	while num > 0:
		hex_str = chars[num % base] + hex_str
		num /= base

	pad_len = width - len(hex_str)
	if pad_len > 0:
		hex_str = '0' * pad_len + hex_str
	
	return hex_str

def from_basex(numstr, chars):
	num = 0
	for i in range(len(numstr)):
		num += len(chars) ** i * chars.index(numstr[-i-1])
	return num

def to_ehex(num, width=0):
	"""Convert number into AIPS ehex."""
	chars = string.digits + string.uppercase
	return to_basex(num, chars, width)

def from_ehex(numstr):
	chars = string.digits + string.uppercase
	return from_basex(numstr, chars)

def get_files(data_dir, id, slot):
	"""Return a list of files that are contained in the dataset."""
	# See Going AIPS C.3 - C.5
	# Based on dircache/regex which proved much faster that globs.
	fn_re = re.compile('(?:CA|CB|GA|HI|MA|PL|SL|UV|AN|BL|BP|CC|CL|CH|CS|FG|FQ|NX|SN|SU).%s...\.%s;' % (to_ehex(slot, width=3), to_ehex(id, width=3)))
	data_files = [data_dir + '/' + file \
		for file in dircache.listdir(data_dir) if fn_re.match(file)]
	return data_files

def data_size(file_list):
	"""Return the total size of the files in bytes."""
	size = 0
	for file in file_list:
		size += os.stat(file).st_size
	return size

def dataset_fs_properties(data_dir, dataset):
	"""Return the owner of the dataset and its size."""
	file_list = get_files(data_dir, dataset['id'], dataset['slot'])
	owner = pwd.getpwuid(os.stat(file_list[0]).st_uid)[0]
	return owner, data_size(file_list)

def num3digit(num):
	"""Convert number to a string with at least 2 signifiant figures."""
	if num < 10:
		return '%.1f' % num
	else:
		return '%3d' % num

def fsize2str(bytes):
	"""Construct a 'pretty' string for an input in bytes."""
	bytes = float(bytes)
	if bytes > 2**40:
		return '%sT' % num3digit(bytes / 2**40)
	elif bytes > 2**30:
		return '%sG' % num3digit(bytes / 2**30)
	elif bytes > 2**20:
		return '%sM' % num3digit(bytes / 2**20)
	elif bytes > 2**10:
		return '%sK' % num3digit(bytes / 2**10)
	else:
		return '%s' % num3digit(bytes)

def str2range(rngstr):
	"""Take a string, eg. '2-5' and return an inclusive range, eg. [2,3,4,5]"""
	subrange = rngstr.split('-')
	try:
		assert len(subrange) == 2
		if subrange[0] == '':
			subrange[0] = '1'
		if subrange[1] == '':
			subrange[1] = str(max_userid)
		assert subrange[0].isdigit() and subrange[1].isdigit()
	except AssertionError:
		raise ValueError, 'parsing "%s"' % rngstr
	return range(int(subrange[0]), int(subrange[1])+1)
	
def parse_range(range_str):
	"""Take a string of values and ranges, eg. '1,2,4-8' and return the numbers."""
	max_userid = 46655	# 36 ** 3 - 1

	range_list = []
	comps = range_str.split(',')
	for rang in comps:
		if rang.isdigit():
			range_list.append(int(rang))
		else:
			try:
				range_list += str2range(rang)
			except ValueError:
				raise ValueError, 'parsing "%s"' % range_str
	return range_list

def aips_disks():
	"""Return list of AIPS disk directories. Disks are
	determined by running DADEVS.SH"""
	data_line = '^\s+Disk .*is (.*)$'
	rx = re.compile(data_line)
	
	try:
		aips_dir = os.environ['AIPS_VERSION']
	except KeyError:
		raise EnvironmentError, 'AIPS_VERSION not defined.'
	cmd = 'sh %s/SYSTEM/UNIX/DADEVS.SH' % aips_dir
	
	status, output = commands.getstatusoutput(cmd)
	if status != 0:
		print >> sys.stderr, 'Problem running DADEVS.SH (%s)' % cmd
		sys.exit(2)
	
	auto_disks = []
	for line in output.split('\n'):
		match = rx.findall(line)
		if match:
			auto_disks.append(match[0])
	if not auto_disks:
		print >> sys.stderr, 'No AIPS data areas found.'
		sys.exit(0)
	
	return auto_disks

def display_dataset(dataset, userids, options, data_dir):
	if dataset['id'] not in userids and 0 not in userids:
		return
	if options.fsize:
		owner, size = dataset_fs_properties(data_dir, dataset)
		fsize_str = fsize2str(size) + ' ' + owner
	else:
		owner, size = None, 0
		fsize_str = ''
	if not options.summary:
		print entry2str(dataset), fsize_str
	return owner, size

def display_disk(data_dir, userids, options):
	if 0 in userids:
		# Display all catalogs
		cat_list = glob.glob('%s/CAD000000.???;' % data_dir)
	else:
		cat_list = []
		for uid in userids:
			cat_file = '%s/CAD000000.%s;' % (data_dir, to_ehex(uid))
			if os.path.exists(cat_file):
				cat_list.append(cat_file)
	cat_list.sort()
	
	if options.headings:
		print data_dir + ':'
		if options.fsize and not options.summary:
			print '  Cat Usid Mapname      Class   Seq  Pt     Last access      Stat Size User'
		elif not options.summary:
			print '  Cat Usid Mapname      Class   Seq  Pt     Last access      Stat'

	user_list = {}
	for cat_file in cat_list:
		catalog = parse_catalog(cat_file)
		for dataset in catalog:
			try:
				owner, size = display_dataset(dataset, userids, options, data_dir)
			except TypeError:
				continue
			if owner in user_list:
				user_list[owner] += size
			else:
				user_list[owner] = size

	if options.headings and options.fsize:
		print 'Total: %s' % fsize2str(sum(user_list.values()))
		for user in user_list:
			print ' %s: %s' % (user, fsize2str(user_list[user]))
		print ''

	return user_list

def main():
	usage = "%prog [options] [dir1] ..."
	parser = optparse.OptionParser(usage=usage, description=__doc__)
	parser.add_option('-a', '--auto-disks', dest='auto_disks', \
			  action='store_true', default=False, \
			  help="""Use AIPS environment to locate data directories.
			  This is the default if no data areas are specified.""")
	parser.add_option('-u', '--userid', metavar='users', \
			  dest='userid', default='0', \
			  help="""Specify user id's. eg. -u 11,14,18-22.
			  Default is to display for all users.""")
	parser.add_option('-t', dest='headings', action='store_false', \
			  default=True, help="Don't display headings and totals.")
	parser.add_option('-l', '--long', dest='fsize', action='store_true', \
			  default=False, help="Display data usage and ownership.")
	parser.add_option('-s', '--summary', dest='summary', action='store_true', \
			  default=False, help="Print disk usage summaries only.")
	(options, disks) = parser.parse_args()
	
	userids = parse_range(options.userid)
	
	if options.summary: 
		options.fsize = True
	if options.summary and not options.headings:
		print >> sys.stderr, 'Can\'t use "summary" and "no titles" modes together.'
		sys.exit(0)
	
	if options.auto_disks or not disks:
		try:
			disks += aips_disks()
		except EnvironmentError:
			if options.auto_disks:
				print >> sys.stderr, '-a specified but $AIPS_VERSION is not defined.'
			else:
				print >> sys.stderr, 'No data areas specified and $AIPS_VERSION not defined.'
			sys.exit(1)
	
	total_user_list = {}
	for data_dir in disks:
		disk_user_list = display_disk(data_dir, userids, options)
		for user in disk_user_list:
			if user in total_user_list:
				total_user_list[user] += disk_user_list[user]
			else:
				total_user_list[user] = disk_user_list[user]

	if options.headings and options.fsize and len(disks) > 1:
		print 'Total (all disks): %s' % fsize2str(sum(total_user_list.values()))
		for user in total_user_list:
			print ' %s: %s' % (user, fsize2str(total_user_list[user]))

if __name__ == '__main__':
	main()
