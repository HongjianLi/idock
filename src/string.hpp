#pragma once
#ifndef IDOCK_STRING_HPP
#define IDOCK_STRING_HPP

#include <string>
#include <istream>
using namespace std;

// Since C++17, copy elision is mandatory and no rvalue reference type or move is required on returning.
// See https://stackoverflow.com/a/38043447 and http://www.open-std.org/jtc1/sc22/wg21/docs/papers/2015/p0135r0.html

//! Remove leading white spaces from a string.
template<class CharT, class Traits, class Allocator>
basic_string<CharT, Traits, Allocator>& trim_start(basic_string<CharT, Traits, Allocator>& str)
{
	decltype(str.size()) i, p = 0;
	for (i = 0; i < str.size() && isspace(str[i]); ++i);
	if (i > 0)
	{
		while (i < str.size()) str[p++] = str[i++];
		str.resize(p);
	}
	return str;
}

//! Remove leading white spaces from a string.
template<class CharT, class Traits, class Allocator>
basic_string<CharT, Traits, Allocator> trim_start(basic_string<CharT, Traits, Allocator>&& str)
{
	decltype(str.size()) i, p = 0;
	for (i = 0; i < str.size() && isspace(str[i]); ++i);
	if (i > 0)
	{
		while (i < str.size()) str[p++] = str[i++];
		str.resize(p);
	}
	return str;
}

//! Remove trailing white spaces from a string.
template<class CharT, class Traits, class Allocator>
basic_string<CharT, Traits, Allocator>& trim_end(basic_string<CharT, Traits, Allocator>& str)
{
	for (auto i = str.size(); i > 0 && isspace(str[i - 1]); --i)
		str.pop_back();
	return str;
}

//! Remove trailing white spaces from a string.
template<class CharT, class Traits, class Allocator>
basic_string<CharT, Traits, Allocator> trim_end(basic_string<CharT, Traits, Allocator>&& str)
{
	for (auto i = str.size(); i > 0 && isspace(str[i - 1]); --i)
		str.pop_back();
	return str;
}

//! Remove leading and trailing white spaces from a string.
template<class CharT, class Traits, class Allocator>
basic_string<CharT, Traits, Allocator>& trim(basic_string<CharT, Traits, Allocator>& str)
{
	return trim_start(trim_end(str));
}

//! Remove leading and trailing white spaces from a string.
template<class CharT, class Traits, class Allocator>
basic_string<CharT, Traits, Allocator> trim(basic_string<CharT, Traits, Allocator>&& str)
{
	return trim_start(trim_end(move(str)));
}

//! Right-aligns the characters by padding them on the left with a specified character, for a specified total length.
template<class CharT, class Traits, class Allocator>
basic_string<CharT, Traits, Allocator>& pad_left(basic_string<CharT, Traits, Allocator>& str, decltype(str.size()) total_width, CharT padding_char)
{
	if (total_width <= str.size())
		return str;
	str.insert(str.begin(), total_width - str.size(), padding_char);
	return str;
}

//! Right-aligns the characters by padding them on the left with a specified character, for a specified total length.
template<class CharT, class Traits, class Allocator>
basic_string<CharT, Traits, Allocator> pad_left(basic_string<CharT, Traits, Allocator>&& str, decltype(str.size()) total_width, CharT padding_char)
{
	if (total_width <= str.size())
		return str;
	str.insert(str.begin(), total_width - str.size(), padding_char);
	return str;
}

//! Left-aligns the characters by padding them on the right with a specified character, for a specified total length.
template<class CharT, class Traits, class Allocator>
basic_string<CharT, Traits, Allocator>& pad_right(basic_string<CharT, Traits, Allocator>& str, decltype(str.size()) total_width, CharT padding_char)
{
	if (total_width <= str.size())
		return str;
	return str.append(total_width - str.size(), padding_char);
}

//! Left-aligns the characters by padding them on the right with a specified character, for a specified total length.
template<class CharT, class Traits, class Allocator>
basic_string<CharT, Traits, Allocator> pad_right(basic_string<CharT, Traits, Allocator>&& str, decltype(str.size()) total_width, CharT padding_char)
{
	if (total_width <= str.size())
		return str;
	return str.append(total_width - str.size(), padding_char);
}

//! Determines whether the beginning of str matches the specified string value.
template<class CharT, class Traits, class Allocator>
bool starts_with(const basic_string<CharT, Traits, Allocator>& str, const basic_string<CharT, Traits, Allocator>& value)
{
	if (str.size() < value.size())
		return false;
	for (auto i1 = value.begin(), i2 = str.begin(); i1 != value.end(); ++i1, ++i2)
		if (*i1 != *i2)
			return false;
	return true;
}

//! Determines whether the beginning of str matches the specified c-style string value.
template<class CharT, class Traits, class Allocator>
bool starts_with(const basic_string<CharT, Traits, Allocator>& str, const CharT* value)
{
	for (auto i = str.begin(); *value; ++i, ++value)
		if (i == str.end() || *i != *value)
			return false;
	return true;
}

//! Determines whether the end of str matches the specified string.
template<class CharT, class Traits, class Allocator>
bool ends_with(const basic_string<CharT, Traits, Allocator>& str, const basic_string<CharT, Traits, Allocator>& value)
{
	if (str.size() < value.size())
		return false;
	for (auto i1 = value.rbegin(), i2 = str.rbegin(); i1 != value.rend(); ++i1, ++i2)
		if (*i1 != *i2)
			return false;
	return true;
}

//! Determines whether the end of str matches the specified c-style string.
template<class CharT, class Traits, class Allocator>
bool ends_with(const basic_string<CharT, Traits, Allocator>& str, const CharT* value)
{
	decltype(str.size()) len = 0;
	for (auto p = value; *p; ++len, ++p);
	if (str.size() < len) return false;
	for (auto i = str.end() - len; *value; ++i, ++value)
		if (i == str.end() || *i != *value)
			return false;
	return true;
}

//! Reads characters from an input stream until any of the Windows (\r\n), Linux (\n) or macOS (\r) line endings is met.
// Origin: https://stackoverflow.com/a/6089413
template<class CharT, class Traits, class Allocator>
basic_istream<CharT, Traits>& safe_getline(basic_istream<CharT, Traits>& input, basic_string<CharT, Traits, Allocator>& str)
{
	str.clear();

	// The characters in the stream are read one-by-one using a std::streambuf.
	// That is faster than reading them one-by-one using the std::istream.
	// Code that uses streambuf this way must be guarded by a sentry object.
	// The sentry object performs various tasks,
	// such as thread synchronization and updating the stream state.

	istream::sentry se(input, true);
	auto sb = input.rdbuf();

	for (;;) {
		int c = sb->sbumpc();
		switch (c) {
		case '\n':
			return input;
		case '\r':
			if (sb->sgetc() == '\n')
				sb->sbumpc();
			return input;
		case streambuf::traits_type::eof():
			// Also handle the case when the last line has no line ending
			if (str.empty())
				input.setstate(ios::eofbit);
			return input;
		default:
			str += (CharT)c;
		}
	}
}

#endif